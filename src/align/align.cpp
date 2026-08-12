#include "align.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <spawn.h>
#include <sstream>
#include <system_error>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#if defined(__linux__)
#include <linux/memfd.h>
#include <sys/syscall.h>
#endif

#include <spdlog/spdlog.h>

extern char** environ;

namespace {

std::string external_insertion_msa_executable;
std::atomic<uint64_t> external_msa_file_counter{0};
std::atomic<uint64_t> external_msa_completed{0};
std::atomic<uint64_t> external_msa_input_nanoseconds{0};
std::atomic<uint64_t> external_msa_process_nanoseconds{0};
std::atomic<uint64_t> external_msa_parse_nanoseconds{0};
std::atomic<uint64_t> external_msa_memfd_inputs{0};
std::atomic<uint64_t> external_msa_file_inputs{0};

constexpr uint32_t kCrossAnchorMinimumInsertion = 10;
constexpr uint32_t kCrossAnchorMaximumDistance = 5;
constexpr uint32_t kCrossAnchorFlankLength = 16;
constexpr double kCrossAnchorMinimumCoverage = 0.70;
constexpr double kCrossAnchorMinimumIdentity = 0.60;
constexpr size_t kCrossAnchorCacheLimit = 1024;

struct CrossAnchorRepairConfiguration {
    std::string executable;
    uint_t maximum_window_span = 3000;
};

CrossAnchorRepairConfiguration cross_anchor_configuration;
std::mutex cross_anchor_configuration_mutex;
std::mutex cross_anchor_cache_mutex;
std::unordered_map<
    std::string,
    std::vector<std::pair<ChrName, std::string>>>
    cross_anchor_msa_cache;

struct CrossAnchorRepairCounters {
    std::atomic<uint64_t> blocks_scanned{0};
    std::atomic<uint64_t> cigar_candidates{0};
    std::atomic<uint64_t> short_insertions_skipped{0};
    std::atomic<uint64_t> same_anchor_skipped{0};
    std::atomic<uint64_t> distance_skipped{0};
    std::atomic<uint64_t> similarity_pairs{0};
    std::atomic<uint64_t> ksw2_repairs{0};
    std::atomic<uint64_t> minipoa_calls{0};
    std::atomic<uint64_t> cache_hits{0};
    std::atomic<uint64_t> accepted{0};
    std::atomic<uint64_t> fallback{0};
    std::atomic<uint64_t> nanoseconds{0};
};

CrossAnchorRepairCounters cross_anchor_counters;

uint64_t elapsedNanoseconds(
    std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point finish) {
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            finish - start)
            .count());
}

std::string ungappedUpper(std::string sequence) {
    sequence.erase(
        std::remove_if(
            sequence.begin(), sequence.end(),
            [](unsigned char c) {
                return c == '-' || std::isspace(c);
            }),
        sequence.end());
    std::transform(
        sequence.begin(), sequence.end(), sequence.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::toupper(c));
        });
    return sequence;
}

class ScopedFd {
public:
    ScopedFd() = default;
    explicit ScopedFd(int fd) : fd_(fd) {}
    ScopedFd(const ScopedFd&) = delete;
    ScopedFd& operator=(const ScopedFd&) = delete;
    ScopedFd(ScopedFd&& other) noexcept : fd_(other.release()) {}
    ScopedFd& operator=(ScopedFd&& other) noexcept {
        if (this != &other) {
            reset(other.release());
        }
        return *this;
    }
    ~ScopedFd() { reset(); }

    int get() const { return fd_; }
    explicit operator bool() const { return fd_ >= 0; }
    int release() {
        const int fd = fd_;
        fd_ = -1;
        return fd;
    }
    void reset(int fd = -1) {
        if (fd_ >= 0) {
            ::close(fd_);
        }
        fd_ = fd;
    }

private:
    int fd_ = -1;
};

class TemporaryMsaInput {
public:
    std::filesystem::path path;

    ~TemporaryMsaInput() {
        std::error_code error;
        if (!path.empty()) {
            std::filesystem::remove(path, error);
        }
    }
};

bool writeAll(int fd, const std::string& contents) {
    size_t offset = 0;
    while (offset < contents.size()) {
        const ssize_t written = ::write(
            fd, contents.data() + offset, contents.size() - offset);
        if (written < 0) {
            if (errno == EINTR) {
                continue;
            }
            return false;
        }
        offset += static_cast<size_t>(written);
    }
    return true;
}

ScopedFd createMsaMemfd(const std::string& contents) {
#if defined(__linux__) && defined(SYS_memfd_create)
    const int fd = static_cast<int>(
        ::syscall(SYS_memfd_create, "ramax-minipoa-input", MFD_CLOEXEC));
    if (fd >= 0) {
        ScopedFd result(fd);
        if (writeAll(fd, contents) && ::lseek(fd, 0, SEEK_SET) == 0) {
            return result;
        }
    }
#else
    (void)contents;
#endif
    return {};
}

bool readAll(int fd, std::string& contents) {
    std::array<char, 8192> buffer{};
    while (true) {
        const ssize_t count = ::read(fd, buffer.data(), buffer.size());
        if (count == 0) {
            return true;
        }
        if (count < 0) {
            if (errno == EINTR) {
                continue;
            }
            return false;
        }
        contents.append(buffer.data(), static_cast<size_t>(count));
    }
}

bool waitForChild(pid_t pid, int& status) {
    while (::waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) {
            return false;
        }
    }
    return true;
}

bool runExternalMsa(
    const std::string& executable,
    std::unordered_map<ChrName, std::string>& sequences) {
    if (executable.empty() ||
        sequences.size() < 2) {
        return false;
    }
    const auto input_start = std::chrono::steady_clock::now();

    std::vector<ChrName> keys;
    keys.reserve(sequences.size());
    for (const auto& [key, unused] : sequences) {
        (void)unused;
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());

    std::ostringstream fasta_stream;
    for (size_t index = 0; index < keys.size(); ++index) {
        fasta_stream << ">s" << index << "\n"
                     << sequences.at(keys[index]) << "\n";
    }
    const std::string fasta = fasta_stream.str();

    ScopedFd input_fd = createMsaMemfd(fasta);
    TemporaryMsaInput temporary;
    std::string input_path;
    if (!input_fd) {
        external_msa_file_inputs.fetch_add(1, std::memory_order_relaxed);
        const uint64_t serial = external_msa_file_counter.fetch_add(
            1, std::memory_order_relaxed);
        const uint64_t stamp = static_cast<uint64_t>(
            std::chrono::steady_clock::now().time_since_epoch().count());
        temporary.path =
            std::filesystem::temp_directory_path() /
            ("ramax-insertion-msa-" + std::to_string(stamp) + "-" +
             std::to_string(serial) + ".input.fa");
        std::ofstream input(temporary.path, std::ios::binary);
        if (!input || !(input << fasta)) {
            throw std::runtime_error(
                "Cannot create external MSA input: " +
                temporary.path.string());
        }
        input_path = temporary.path.string();
    } else {
        external_msa_memfd_inputs.fetch_add(1, std::memory_order_relaxed);
    }

    int output_pipe[2] = {-1, -1};
#if defined(__linux__)
    const int pipe_result = ::pipe2(output_pipe, O_CLOEXEC);
#else
    const int pipe_result = ::pipe(output_pipe);
#endif
    if (pipe_result != 0) {
        throw std::system_error(
            errno, std::generic_category(), "Cannot create MSA output pipe");
    }
    ScopedFd output_read(output_pipe[0]);
    ScopedFd output_write(output_pipe[1]);
#if !defined(__linux__)
    for (const int fd : output_pipe) {
        const int flags = ::fcntl(fd, F_GETFD);
        if (flags < 0 || ::fcntl(fd, F_SETFD, flags | FD_CLOEXEC) < 0) {
            throw std::system_error(
                errno, std::generic_category(),
                "Cannot mark MSA output pipe close-on-exec");
        }
    }
#endif

    int child_input_fd = 100;
    while (child_input_fd == input_fd.get() ||
           child_input_fd == output_read.get() ||
           child_input_fd == output_write.get()) {
        ++child_input_fd;
    }
    if (input_fd) {
        input_path =
            "/proc/self/fd/" + std::to_string(child_input_fd);
    }

    posix_spawn_file_actions_t actions;
    if (posix_spawn_file_actions_init(&actions) != 0) {
        throw std::runtime_error("Cannot initialize MSA spawn actions");
    }
    const auto destroy_actions = [&]() {
        posix_spawn_file_actions_destroy(&actions);
    };
    int action_error = 0;
    if (input_fd) {
        action_error = posix_spawn_file_actions_adddup2(
            &actions, input_fd.get(), child_input_fd);
    }
    if (action_error == 0) {
        action_error = posix_spawn_file_actions_adddup2(
            &actions, output_write.get(), STDOUT_FILENO);
    }
    if (action_error == 0) {
        action_error = posix_spawn_file_actions_addclose(
            &actions, output_read.get());
    }
    if (action_error == 0 && output_write.get() != STDOUT_FILENO) {
        action_error = posix_spawn_file_actions_addclose(
            &actions, output_write.get());
    }
    if (action_error != 0) {
        destroy_actions();
        throw std::system_error(
            action_error, std::generic_category(),
            "Cannot configure MSA spawn actions");
    }

    std::vector<char*> arguments{
        const_cast<char*>(executable.c_str()),
        const_cast<char*>("-r"),
        const_cast<char*>("1"),
        const_cast<char*>("-t"),
        const_cast<char*>("1"),
        input_path.data(),
        nullptr};
    pid_t pid = -1;
    const auto process_start = std::chrono::steady_clock::now();
    external_msa_input_nanoseconds.fetch_add(
        elapsedNanoseconds(input_start, process_start),
        std::memory_order_relaxed);
    const int spawn_error = posix_spawnp(
        &pid, executable.c_str(), &actions, nullptr,
        arguments.data(), environ);
    destroy_actions();
    if (spawn_error != 0) {
        spdlog::warn(
            "[external-msa] cannot spawn {}: {}",
            executable, std::generic_category().message(spawn_error));
        return false;
    }
    output_write.reset();

    std::string output_text;
    const bool output_read_ok = readAll(output_read.get(), output_text);
    output_read.reset();
    int child_status = 0;
    const bool waited = waitForChild(pid, child_status);
    if (!output_read_ok || !waited || !WIFEXITED(child_status) ||
        WEXITSTATUS(child_status) != 0) {
        spdlog::warn(
            "[external-msa] command failed: {}",
            executable);
        return false;
    }

    const auto parse_start = std::chrono::steady_clock::now();
    external_msa_process_nanoseconds.fetch_add(
        elapsedNanoseconds(process_start, parse_start),
        std::memory_order_relaxed);
    std::istringstream output(output_text);

    std::unordered_map<std::string, std::string> aligned_by_id;
    std::string current_id;
    std::string line;
    while (std::getline(output, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line.front() == '>') {
            std::istringstream header(line.substr(1));
            header >> current_id;
            if (current_id.empty() ||
                !aligned_by_id.emplace(current_id, std::string{}).second) {
                spdlog::warn(
                    "[external-msa] invalid FASTA header");
                return false;
            }
            continue;
        }
        if (current_id.empty()) {
            if (line.empty()) {
                continue;
            }
            spdlog::warn(
                "[external-msa] sequence before FASTA header");
            return false;
        }
        for (const unsigned char c : line) {
            if (!std::isspace(c)) {
                aligned_by_id[current_id].push_back(
                    static_cast<char>(c));
            }
        }
    }

    if (aligned_by_id.size() != keys.size()) {
        spdlog::warn("[external-msa] output row count mismatch");
        return false;
    }

    size_t aligned_length = 0;
    std::vector<std::string> aligned(keys.size());
    for (size_t index = 0; index < keys.size(); ++index) {
        const std::string id = "s" + std::to_string(index);
        const auto aligned_it = aligned_by_id.find(id);
        if (aligned_it == aligned_by_id.end() ||
            ungappedUpper(aligned_it->second) !=
                ungappedUpper(sequences.at(keys[index]))) {
            spdlog::warn(
                "[external-msa] output sequence validation failed");
            return false;
        }
        if (index == 0) {
            aligned_length = aligned_it->second.size();
        } else if (aligned_it->second.size() != aligned_length) {
            spdlog::warn(
                "[external-msa] output rows have unequal lengths");
            return false;
        }
        aligned[index] = aligned_it->second;
    }

    for (size_t index = 0; index < keys.size(); ++index) {
        sequences[keys[index]] = std::move(aligned[index]);
    }

    const auto parse_finish = std::chrono::steady_clock::now();
    external_msa_parse_nanoseconds.fetch_add(
        elapsedNanoseconds(parse_start, parse_finish),
        std::memory_order_relaxed);

    const uint64_t completed =
        external_msa_completed.fetch_add(1, std::memory_order_relaxed) + 1;
    if (completed % 1000 == 0) {
        constexpr double kNanosecondsPerSecond = 1.0e9;
        spdlog::info(
            "[external-msa] completed={} executable={} input_seconds={:.3f} "
            "process_seconds={:.3f} parse_seconds={:.3f} memfd_inputs={} "
            "file_inputs={}",
            completed, executable,
            external_msa_input_nanoseconds.load(
                std::memory_order_relaxed) /
                kNanosecondsPerSecond,
            external_msa_process_nanoseconds.load(
                std::memory_order_relaxed) /
                kNanosecondsPerSecond,
            external_msa_parse_nanoseconds.load(
                std::memory_order_relaxed) /
                kNanosecondsPerSecond,
            external_msa_memfd_inputs.load(std::memory_order_relaxed),
            external_msa_file_inputs.load(std::memory_order_relaxed));
    }
    return true;
}

}  // namespace

static uint_t mergeAlignmentByRefCore(
    const ChrName& ref_name,
    std::unordered_map<ChrName, std::string>& seqs,
    const std::unordered_map<ChrName, Cigar_t>& cigars);

bool alignSequencesWithExternalMsa(
    const std::string& executable,
    std::unordered_map<ChrName, std::string>& sequences) {
    return runExternalMsa(executable, sequences);
}

void configureExternalInsertionMsa(const std::string& executable) {
    external_insertion_msa_executable = executable;
    if (!executable.empty()) {
        spdlog::info(
            "[external-msa] insertion aligner configured: {}",
            executable);
    }
}

void InsertInfo::alignSeqs() {
    if (aligned || seqs.empty()) {
        return;
    }

    if (seqs.size() == 1) {
        ref_name = seqs.begin()->first;
        total_length = seqs.begin()->second.size();
        aligned = true;
        return;
    }

    if (runExternalMsa(
            external_insertion_msa_executable, seqs)) {
        ref_name = seqs.begin()->first;
        total_length = seqs.begin()->second.size();
        aligned = true;
        return;
    }

    size_t max_len = 0;
    ChrName longest_key;
    for (const auto& [key, sequence] : seqs) {
        if (sequence.size() > max_len) {
            max_len = sequence.size();
            longest_key = key;
        }
    }
    ref_name = longest_key;

    std::unordered_map<ChrName, Cigar_t> cigars;
    for (const auto& [key, sequence] : seqs) {
        (void)sequence;
        if (key == ref_name) {
            continue;
        }
        cigars[key] = globalAlignKSW2(seqs[ref_name], seqs.at(key));
    }
    total_length = mergeAlignmentByRefCore(ref_name, seqs, cigars);
    aligned = true;
}

KSW2AlignConfig makeDefaultKSW2Config() {
    static int8_t simple_dna_mat[25];
    static bool initialized = false;
    if (!initialized) {
        // A C G T N -> 0 1 2 3 4
        const int match = 2;
        const int mismatch = -3;
        const int ambiguous = -1;  // 对N的惩罚较小

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 4 || j == 4) {
                    simple_dna_mat[i * 5 + j] = ambiguous;  // N匹配
                }
                else if (i == j) {
                    simple_dna_mat[i * 5 + j] = match;
                }
                else {
                    simple_dna_mat[i * 5 + j] = mismatch;
                }
            }
        }
        initialized = true;
    }

    return {
        .mat = simple_dna_mat,
        .alphabet_size = 5,
        .gap_open = 8,
        .gap_extend = 1,
        .end_bonus = 0,
        .zdrop = 100,           // 较远匹配终止
        .band_width = -1,      // 合理的band可提高性能
        .flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT
    };
}

Cigar_t globalAlignKSW2(const std::string& ref,
    const std::string& query)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    init_simd_mat();

    KSW2AlignConfig cfg;
    cfg.mat = dna5_simd_mat;
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;          // gap open penalty
    cfg.gap_extend = 2;        // gap extension penalty
    cfg.end_bonus = 0;         // ❌ 不需要 ends-free 奖励
    cfg.zdrop = -1;            // ❌ 禁用 z-drop（全局比对必须完整比完）
    cfg.band_width = -1;       // 启用全矩阵（也可设 auto_band）
    cfg.flag = KSW_EZ_RIGHT; // ✅ 通用矩阵 + gap右对齐

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}

Cigar_t globalAlignKSW2_2(const std::string& ref,
    const std::string& query)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    KSW2AlignConfig cfg = makeTurboKSW2Config2(query.size(), ref.size());

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}


/**********************************************************************
*  extendAlignKSW2  ——  ends-free（seed-and-extend）比对
*    @param ref        参考片段（目标方向）
*    @param query      查询片段（同方向；若反链请先反向互补）
*    @param zdrop      Z-drop 剪枝阈值（默认 200）
*    @param band       带宽限制；<0 表示不限制
*    @return           Cigar_t（BAM 编码）
**********************************************************************/
Cigar_t extendAlignKSW2(const std::string& ref,
    const std::string& query,
    int zdrop)
{
    /* ---------- 1. 序列编码 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());
    for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = ScoreChar2Idx[(uint8_t)ref[i]];
    for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = ScoreChar2Idx[(uint8_t)query[i]];

    ///* ---------- 2. 配置 ---------- */
    //KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
    ////KSW2AlignConfig cfg;
    //cfg.zdrop = zdrop;       // 用于提前终止
    //cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
    //    | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
    //    | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
    //    | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
    //// 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
    //cfg.end_bonus = 100;
    //cfg.band_width = -1;
    init_simd_mat();
    KSW2AlignConfig cfg;
	cfg.mat = dna5_simd_mat;
    cfg.zdrop = zdrop;
    cfg.flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_APPROX_DROP;
    cfg.end_bonus = 50;
    cfg.alphabet_size = 5;
    cfg.gap_open = 5;
    cfg.gap_extend = 2;
    cfg.band_width = auto_band(ref.size(), query.size());


    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};
    ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);

    // 赋值bool& if_zdrop,int& ref_end,int& qry_end
    /* ---------- 4. 拷贝 & 释放 ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);                    // ksw2 使用 malloc
    return cigar;                      // 返回的 CIGAR 即延伸片段
}

/* ──────────── 合并成 MSA (就地修改 seqs) ──────────── */
static uint_t mergeAlignmentByRefCore(
    const ChrName& ref_name,
    std::unordered_map<ChrName, std::string>& seqs,
    const std::unordered_map<ChrName, Cigar_t>& cigars)
{
    auto ref_it = seqs.find(ref_name);
    if (ref_it == seqs.end())
        throw std::invalid_argument("mergeAlignmentByRef: ref not found");

    std::string& ref_raw = ref_it->second;
    uint_t total_aligned_length = ref_raw.size();
	RefAlignInfo insert_info;

    for (const auto& [key, cigar] : cigars) {
        if (key == ref_name) continue;
        auto q_it = seqs.find(key);
        if (q_it == seqs.end()) {
            throw std::invalid_argument("mergeAlignmentByRef: seq missing");
        }

		std::string& qry_raw = q_it->second;

		uint_t ref_pos = 0;
		uint_t qry_pos = 0;
        for (auto& unit : cigar) {
            uint32_t len;
            char op;
			intToCigar(unit, op, len);

            if (op == 'D') {
				qry_raw.insert(qry_pos, len, '-');
                ref_pos += len;
                qry_pos += len;
            }
            else if (op == 'I') {
                std::string ins = qry_raw.substr(qry_pos, len);

                // 2) 在 insert_info 里插入或更新
                auto it = insert_info.find(ref_pos);
                if (it != insert_info.end()) {
                    it->second.seqs[key] = ins;
                }
                else {
                    InsertInfo info;
                    info.seqs[key] = ins;
                    insert_info[ref_pos] = std::move(info);
                }

                // 3) 从原始 query 序列里移除这段已“消费”的子串
                qry_raw.erase(qry_pos, len);
            }
            else {
                ref_pos += len;
                qry_pos += len;
            }
        }
    }

    uint_t offset = 0;
	for (auto& [ref_pos, info] : insert_info) {
		info.alignSeqs(); // 对齐所有插入序列
		if (info.ref_name.empty()) continue; // 没有参考序列，跳过
		for (auto& [sp_name, seq] : seqs) {
			auto it = info.seqs.find(sp_name);
            if (it != info.seqs.end()) {
				seq.insert(ref_pos + offset, it->second); // 在 ref_pos 位置插入
            }
            else {
                seq.insert(ref_pos + offset, info.total_length, '-');   // 直接用 string::insert 重载
            }
		}
        offset += info.total_length; // 更新总长度
        total_aligned_length += info.total_length;
	}

    for (auto& [chr, seq] : seqs) {
        if (seq.size() != total_aligned_length) {
            std::cout << "";
        }
    }

    return total_aligned_length;

}

namespace {

struct CrossAnchorInsertionEvent {
    ChrName key;
    uint32_t reference_position = 0;
    uint32_t query_position = 0;
    uint32_t length = 0;
};

struct CrossAnchorSimilarity {
    size_t left = 0;
    size_t right = 0;
    size_t paired_bases = 0;
    size_t matching_bases = 0;
    double coverage = 0.0;
    double identity = 0.0;
};

struct CrossAnchorGroup {
    std::vector<size_t> event_indices;
    std::vector<CrossAnchorSimilarity> edges;
    uint32_t minimum_position = 0;
    uint32_t maximum_position = 0;
    uint32_t left_reference = 0;
    uint32_t right_reference = 0;
    size_t species_count = 0;
    size_t shared_bases = 0;
    double average_identity = 0.0;
};

struct PairAlignmentMetrics {
    size_t paired_bases = 0;
    size_t matching_bases = 0;
};

char normalizedBase(char base) {
    return static_cast<char>(
        std::toupper(static_cast<unsigned char>(base)));
}

std::string removeAlignmentGaps(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (const char base : sequence) {
        if (base != '-') {
            result.push_back(base);
        }
    }
    return result;
}

Cigar_t globalAlignKSW2Banded(const std::string& reference,
                              const std::string& query) {
    if (reference.empty()) {
        Cigar_t result;
        if (!query.empty()) {
            appendCigarOp(
                result, 'I', static_cast<uint32_t>(query.size()));
        }
        return result;
    }
    if (query.empty()) {
        Cigar_t result;
        appendCigarOp(
            result, 'D', static_cast<uint32_t>(reference.size()));
        return result;
    }

    std::vector<uint8_t> reference_encoded(reference.size());
    std::vector<uint8_t> query_encoded(query.size());
    for (size_t index = 0; index < reference.size(); ++index) {
        reference_encoded[index] =
            ScoreChar2Idx[static_cast<uint8_t>(reference[index])];
    }
    for (size_t index = 0; index < query.size(); ++index) {
        query_encoded[index] =
            ScoreChar2Idx[static_cast<uint8_t>(query[index])];
    }

    init_simd_mat();
    const size_t difference = reference.size() > query.size()
        ? reference.size() - query.size()
        : query.size() - reference.size();
    const size_t requested_band = difference + 32;
    const int band = requested_band >
            static_cast<size_t>(std::numeric_limits<int>::max())
        ? std::numeric_limits<int>::max()
        : static_cast<int>(requested_band);

    ksw_extz_t alignment{};
    ksw_extz2_sse(
        nullptr,
        static_cast<int>(query_encoded.size()), query_encoded.data(),
        static_cast<int>(reference_encoded.size()), reference_encoded.data(),
        5, dna5_simd_mat, 5, 2, band, -1, 0,
        KSW_EZ_RIGHT, &alignment);

    Cigar_t result;
    result.reserve(alignment.n_cigar);
    for (int index = 0; index < alignment.n_cigar; ++index) {
        result.push_back(alignment.cigar[index]);
    }
    free(alignment.cigar);
    return result;
}

std::optional<CrossAnchorSimilarity> insertionSimilarity(
    size_t left_index,
    size_t right_index,
    const CrossAnchorInsertionEvent& left_event,
    const CrossAnchorInsertionEvent& right_event,
    const std::unordered_map<ChrName, std::string>& sequences) {
    const auto left_sequence_it = sequences.find(left_event.key);
    const auto right_sequence_it = sequences.find(right_event.key);
    if (left_sequence_it == sequences.end() ||
        right_sequence_it == sequences.end() ||
        left_event.query_position + left_event.length >
            left_sequence_it->second.size() ||
        right_event.query_position + right_event.length >
            right_sequence_it->second.size()) {
        return std::nullopt;
    }

    const std::string left_sequence = left_sequence_it->second.substr(
        left_event.query_position, left_event.length);
    const std::string right_sequence = right_sequence_it->second.substr(
        right_event.query_position, right_event.length);
    const Cigar_t cigar = globalAlignKSW2Banded(
        left_sequence, right_sequence);

    size_t left_position = 0;
    size_t right_position = 0;
    size_t paired = 0;
    size_t matching = 0;
    for (const auto unit : cigar) {
        uint32_t length = 0;
        char operation = '\0';
        intToCigar(unit, operation, length);
        if (operation == 'M' || operation == '=' || operation == 'X') {
            if (left_position + length > left_sequence.size() ||
                right_position + length > right_sequence.size()) {
                return std::nullopt;
            }
            for (uint32_t offset = 0; offset < length; ++offset) {
                matching += normalizedBase(
                    left_sequence[left_position + offset]) ==
                    normalizedBase(
                        right_sequence[right_position + offset]);
            }
            paired += length;
            left_position += length;
            right_position += length;
        } else if (operation == 'I') {
            right_position += length;
        } else if (operation == 'D') {
            left_position += length;
        } else {
            return std::nullopt;
        }
    }
    if (left_position != left_sequence.size() ||
        right_position != right_sequence.size() || paired == 0) {
        return std::nullopt;
    }

    const size_t shorter = std::min(
        left_sequence.size(), right_sequence.size());
    const double coverage = shorter == 0
        ? 0.0
        : static_cast<double>(paired) /
            static_cast<double>(shorter);
    const double identity = static_cast<double>(matching) /
        static_cast<double>(paired);
    if (paired < kCrossAnchorMinimumInsertion ||
        coverage < kCrossAnchorMinimumCoverage ||
        identity < kCrossAnchorMinimumIdentity) {
        return std::nullopt;
    }
    return CrossAnchorSimilarity{
        left_index, right_index, paired, matching, coverage, identity};
}

PairAlignmentMetrics alignedPairMetrics(const std::string& left,
                                        const std::string& right) {
    PairAlignmentMetrics result;
    if (left.size() != right.size()) {
        return result;
    }
    for (size_t column = 0; column < left.size(); ++column) {
        if (left[column] == '-' || right[column] == '-') {
            continue;
        }
        ++result.paired_bases;
        result.matching_bases +=
            normalizedBase(left[column]) == normalizedBase(right[column]);
    }
    return result;
}

std::vector<CrossAnchorInsertionEvent> collectInsertionEvents(
    const ChrName& reference_name,
    const std::unordered_map<ChrName, std::string>& sequences,
    const std::unordered_map<ChrName, Cigar_t>& cigars) {
    std::vector<CrossAnchorInsertionEvent> events;
    for (const auto& [key, cigar] : cigars) {
        if (key == reference_name || sequences.find(key) == sequences.end()) {
            continue;
        }
        uint64_t reference_position = 0;
        uint64_t query_position = 0;
        for (const auto unit : cigar) {
            uint32_t length = 0;
            char operation = '\0';
            intToCigar(unit, operation, length);
            if (operation == 'M' || operation == '=' || operation == 'X') {
                reference_position += length;
                query_position += length;
            } else if (operation == 'I') {
                if (length >= kCrossAnchorMinimumInsertion &&
                    reference_position <=
                        std::numeric_limits<uint32_t>::max() &&
                    query_position <=
                        std::numeric_limits<uint32_t>::max()) {
                    events.push_back(CrossAnchorInsertionEvent{
                        key,
                        static_cast<uint32_t>(reference_position),
                        static_cast<uint32_t>(query_position),
                        length});
                } else if (length < kCrossAnchorMinimumInsertion) {
                    cross_anchor_counters.short_insertions_skipped.fetch_add(
                        1, std::memory_order_relaxed);
                }
                query_position += length;
            } else if (operation == 'D' || operation == 'N') {
                reference_position += length;
            } else if (operation == 'S') {
                query_position += length;
            }
        }
    }
    return events;
}

class ConstrainedDisjointSet {
public:
    explicit ConstrainedDisjointSet(
        const std::vector<CrossAnchorInsertionEvent>& events)
        : parent_(events.size()), rank_(events.size(), 0),
          minimum_(events.size()), maximum_(events.size()) {
        std::iota(parent_.begin(), parent_.end(), 0);
        for (size_t index = 0; index < events.size(); ++index) {
            minimum_[index] = events[index].reference_position;
            maximum_[index] = events[index].reference_position;
        }
    }

    size_t find(size_t value) {
        if (parent_[value] != value) {
            parent_[value] = find(parent_[value]);
        }
        return parent_[value];
    }

    bool unite(size_t left, size_t right) {
        left = find(left);
        right = find(right);
        if (left == right) {
            return true;
        }
        const uint32_t merged_minimum =
            std::min(minimum_[left], minimum_[right]);
        const uint32_t merged_maximum =
            std::max(maximum_[left], maximum_[right]);
        if (merged_maximum - merged_minimum >
            kCrossAnchorMaximumDistance) {
            return false;
        }
        if (rank_[left] < rank_[right]) {
            std::swap(left, right);
        }
        parent_[right] = left;
        minimum_[left] = merged_minimum;
        maximum_[left] = merged_maximum;
        if (rank_[left] == rank_[right]) {
            ++rank_[left];
        }
        return true;
    }

private:
    std::vector<size_t> parent_;
    std::vector<uint8_t> rank_;
    std::vector<uint32_t> minimum_;
    std::vector<uint32_t> maximum_;
};

std::vector<CrossAnchorGroup> buildCrossAnchorGroups(
    const std::vector<CrossAnchorInsertionEvent>& events,
    const std::unordered_map<ChrName, std::string>& sequences) {
    if (events.size() < 2) {
        return {};
    }

    std::vector<CrossAnchorSimilarity> edges;
    bool structural_candidate = false;
    for (size_t left = 0; left < events.size(); ++left) {
        for (size_t right = left + 1; right < events.size(); ++right) {
            if (events[left].key == events[right].key) {
                continue;
            }
            const uint32_t left_position = events[left].reference_position;
            const uint32_t right_position = events[right].reference_position;
            if (left_position == right_position) {
                cross_anchor_counters.same_anchor_skipped.fetch_add(
                    1, std::memory_order_relaxed);
                continue;
            }
            const uint32_t distance = left_position > right_position
                ? left_position - right_position
                : right_position - left_position;
            if (distance > kCrossAnchorMaximumDistance) {
                cross_anchor_counters.distance_skipped.fetch_add(
                    1, std::memory_order_relaxed);
                continue;
            }
            structural_candidate = true;
            const auto similarity = insertionSimilarity(
                left, right, events[left], events[right], sequences);
            if (similarity) {
                edges.push_back(*similarity);
                cross_anchor_counters.similarity_pairs.fetch_add(
                    1, std::memory_order_relaxed);
            }
        }
    }
    if (structural_candidate) {
        cross_anchor_counters.cigar_candidates.fetch_add(
            1, std::memory_order_relaxed);
    }
    if (edges.empty()) {
        return {};
    }

    std::sort(edges.begin(), edges.end(),
        [&](const auto& left, const auto& right) {
            if (left.paired_bases != right.paired_bases) {
                return left.paired_bases > right.paired_bases;
            }
            if (left.identity != right.identity) {
                return left.identity > right.identity;
            }
            const uint32_t left_coordinate = std::min(
                events[left.left].reference_position,
                events[left.right].reference_position);
            const uint32_t right_coordinate = std::min(
                events[right.left].reference_position,
                events[right.right].reference_position);
            return left_coordinate < right_coordinate;
        });

    ConstrainedDisjointSet components(events);
    std::vector<CrossAnchorSimilarity> accepted_edges;
    for (const auto& edge : edges) {
        if (components.unite(edge.left, edge.right)) {
            accepted_edges.push_back(edge);
        }
    }

    std::map<size_t, CrossAnchorGroup> grouped;
    for (size_t index = 0; index < events.size(); ++index) {
        const size_t root = components.find(index);
        grouped[root].event_indices.push_back(index);
    }
    for (const auto& edge : accepted_edges) {
        const size_t root = components.find(edge.left);
        if (components.find(edge.right) == root) {
            grouped[root].edges.push_back(edge);
        }
    }

    std::vector<CrossAnchorGroup> result;
    for (auto& [root, group] : grouped) {
        (void)root;
        if (group.edges.empty()) {
            continue;
        }
        std::set<ChrName> species;
        group.minimum_position = std::numeric_limits<uint32_t>::max();
        group.maximum_position = 0;
        double identity_sum = 0.0;
        for (const size_t index : group.event_indices) {
            species.insert(events[index].key);
            group.minimum_position = std::min(
                group.minimum_position,
                events[index].reference_position);
            group.maximum_position = std::max(
                group.maximum_position,
                events[index].reference_position);
        }
        if (species.size() < 2 ||
            group.minimum_position == group.maximum_position ||
            group.maximum_position - group.minimum_position >
                kCrossAnchorMaximumDistance) {
            continue;
        }
        for (const auto& edge : group.edges) {
            group.shared_bases += edge.paired_bases;
            identity_sum += edge.identity;
        }
        group.species_count = species.size();
        group.average_identity = identity_sum /
            static_cast<double>(group.edges.size());
        result.push_back(std::move(group));
    }
    return result;
}

std::string buildCrossAnchorCacheKey(
    const std::string& executable,
    const std::map<ChrName, std::string>& alleles) {
    std::string key;
    key.reserve(executable.size() + 64);
    key.append(executable);
    key.push_back('\0');
    for (const auto& [name, sequence] : alleles) {
        key.append(std::to_string(name.size()));
        key.push_back(':');
        key.append(name);
        key.append(std::to_string(sequence.size()));
        key.push_back(':');
        key.append(sequence);
        key.push_back('\0');
    }
    return key;
}

bool alignLocalAlleles(
    const std::map<ChrName, std::string>& raw_alleles,
    size_t insertion_species_count,
    const CrossAnchorRepairConfiguration& configuration,
    std::unordered_map<ChrName, std::string>& aligned) {
    if (raw_alleles.size() < 2) {
        return false;
    }
    aligned.clear();
    for (const auto& [key, sequence] : raw_alleles) {
        aligned.emplace(key, sequence);
    }

    if (insertion_species_count >= 3) {
        if (configuration.executable.empty()) {
            return false;
        }
        const std::string cache_key = buildCrossAnchorCacheKey(
            configuration.executable, raw_alleles);
        {
            std::lock_guard lock(cross_anchor_cache_mutex);
            const auto cached = cross_anchor_msa_cache.find(cache_key);
            if (cached != cross_anchor_msa_cache.end()) {
                aligned.clear();
                for (const auto& [key, sequence] : cached->second) {
                    aligned.emplace(key, sequence);
                }
                cross_anchor_counters.cache_hits.fetch_add(
                    1, std::memory_order_relaxed);
                return true;
            }
        }
        cross_anchor_counters.minipoa_calls.fetch_add(
            1, std::memory_order_relaxed);
        if (!runExternalMsa(configuration.executable, aligned)) {
            return false;
        }
        std::vector<std::pair<ChrName, std::string>> cached_rows;
        cached_rows.reserve(aligned.size());
        for (const auto& [key, sequence] : aligned) {
            cached_rows.emplace_back(key, sequence);
        }
        std::sort(cached_rows.begin(), cached_rows.end());
        {
            std::lock_guard lock(cross_anchor_cache_mutex);
            if (cross_anchor_msa_cache.size() >= kCrossAnchorCacheLimit) {
                cross_anchor_msa_cache.clear();
            }
            cross_anchor_msa_cache.emplace(
                cache_key, std::move(cached_rows));
        }
        return true;
    }

    auto longest = raw_alleles.begin();
    for (auto iterator = std::next(raw_alleles.begin());
         iterator != raw_alleles.end(); ++iterator) {
        if (iterator->second.size() > longest->second.size()) {
            longest = iterator;
        }
    }
    std::unordered_map<ChrName, Cigar_t> local_cigars;
    for (const auto& [key, sequence] : raw_alleles) {
        if (key != longest->first) {
            local_cigars.emplace(
                key,
                globalAlignKSW2Banded(longest->second, sequence));
        }
    }
    try {
        mergeAlignmentByRefCore(longest->first, aligned, local_cigars);
    } catch (const std::exception&) {
        return false;
    }
    cross_anchor_counters.ksw2_repairs.fetch_add(
        1, std::memory_order_relaxed);
    return true;
}

bool rangesOverlap(const CrossAnchorGroup& left,
                   const CrossAnchorGroup& right) {
    return left.left_reference < right.right_reference &&
        right.left_reference < left.right_reference;
}

struct ReferenceAlignmentLayout {
    std::vector<std::pair<size_t, size_t>> gap_ranges;
    std::vector<size_t> base_columns;
};

std::optional<ReferenceAlignmentLayout> buildReferenceLayout(
    const std::string& aligned_reference) {
    ReferenceAlignmentLayout layout;
    size_t gap_start = 0;
    for (size_t column = 0; column < aligned_reference.size(); ++column) {
        if (aligned_reference[column] == '-') {
            continue;
        }
        layout.gap_ranges.emplace_back(gap_start, column);
        layout.base_columns.push_back(column);
        gap_start = column + 1;
    }
    layout.gap_ranges.emplace_back(gap_start, aligned_reference.size());
    if (layout.gap_ranges.size() != layout.base_columns.size() + 1) {
        return std::nullopt;
    }
    return layout;
}

bool mergeRepairedAndLegacyProfiles(
    const ChrName& reference_name,
    const std::set<ChrName>& participant_keys,
    const std::unordered_map<ChrName, std::string>& repaired,
    const std::unordered_map<ChrName, std::string>& legacy,
    std::unordered_map<ChrName, std::string>& merged) {
    const auto repaired_reference_it = repaired.find(reference_name);
    const auto legacy_reference_it = legacy.find(reference_name);
    if (repaired_reference_it == repaired.end() ||
        legacy_reference_it == legacy.end()) {
        return false;
    }
    const auto repaired_layout = buildReferenceLayout(
        repaired_reference_it->second);
    const auto legacy_layout = buildReferenceLayout(
        legacy_reference_it->second);
    if (!repaired_layout || !legacy_layout ||
        repaired_layout->base_columns.size() !=
            legacy_layout->base_columns.size()) {
        return false;
    }

    for (size_t index = 0;
         index < repaired_layout->base_columns.size(); ++index) {
        const char repaired_base = repaired_reference_it->second[
            repaired_layout->base_columns[index]];
        const char legacy_base = legacy_reference_it->second[
            legacy_layout->base_columns[index]];
        if (normalizedBase(repaired_base) != normalizedBase(legacy_base)) {
            return false;
        }
    }

    merged.clear();
    for (const auto& [key, sequence] : legacy) {
        (void)sequence;
        merged.emplace(key, std::string{});
    }
    for (const auto& [key, sequence] : repaired) {
        (void)sequence;
        if (merged.find(key) == merged.end()) {
            return false;
        }
    }

    const size_t reference_bases = repaired_layout->base_columns.size();
    for (size_t boundary = 0; boundary <= reference_bases; ++boundary) {
        std::vector<size_t> repaired_gap_columns;
        const auto [repaired_begin, repaired_end] =
            repaired_layout->gap_ranges[boundary];
        for (size_t column = repaired_begin;
             column < repaired_end; ++column) {
            const bool occupied = std::any_of(
                participant_keys.begin(), participant_keys.end(),
                [&](const ChrName& key) {
                    const auto row = repaired.find(key);
                    return row != repaired.end() &&
                        column < row->second.size() &&
                        row->second[column] != '-';
                });
            if (occupied) {
                repaired_gap_columns.push_back(column);
            }
        }

        std::vector<size_t> legacy_gap_columns;
        const auto [legacy_begin, legacy_end] =
            legacy_layout->gap_ranges[boundary];
        for (size_t column = legacy_begin;
             column < legacy_end; ++column) {
            bool occupied = false;
            for (const auto& [key, row] : legacy) {
                if (key != reference_name &&
                    participant_keys.find(key) == participant_keys.end() &&
                    column < row.size() && row[column] != '-') {
                    occupied = true;
                    break;
                }
            }
            if (occupied) {
                legacy_gap_columns.push_back(column);
            }
        }

        for (const size_t column : repaired_gap_columns) {
            for (auto& [key, row] : merged) {
                if (key == reference_name ||
                    participant_keys.find(key) == participant_keys.end()) {
                    row.push_back('-');
                } else {
                    row.push_back(repaired.at(key)[column]);
                }
            }
        }
        for (const size_t column : legacy_gap_columns) {
            for (auto& [key, row] : merged) {
                if (key == reference_name ||
                    participant_keys.find(key) != participant_keys.end()) {
                    row.push_back('-');
                } else {
                    row.push_back(legacy.at(key)[column]);
                }
            }
        }

        if (boundary == reference_bases) {
            continue;
        }
        const size_t repaired_column =
            repaired_layout->base_columns[boundary];
        const size_t legacy_column = legacy_layout->base_columns[boundary];
        for (auto& [key, row] : merged) {
            if (key == reference_name ||
                participant_keys.find(key) != participant_keys.end()) {
                row.push_back(repaired.at(key)[repaired_column]);
            } else {
                row.push_back(legacy.at(key)[legacy_column]);
            }
        }
    }
    return true;
}

bool repairCrossAnchorGroups(
    const ChrName& reference_name,
    size_t raw_reference_length,
    const std::vector<CrossAnchorInsertionEvent>& events,
    std::vector<CrossAnchorGroup> groups,
    const CrossAnchorRepairConfiguration& configuration,
    std::unordered_map<ChrName, std::string>& sequences) {
    const auto reference_it = sequences.find(reference_name);
    if (reference_it == sequences.end()) {
        return false;
    }
    const std::string& aligned_reference = reference_it->second;
    std::vector<size_t> reference_boundaries(
        raw_reference_length + 1, aligned_reference.size());
    reference_boundaries[0] = 0;
    size_t reference_position = 0;
    for (size_t column = 0; column < aligned_reference.size(); ++column) {
        if (aligned_reference[column] != '-') {
            if (++reference_position > raw_reference_length) {
                return false;
            }
            reference_boundaries[reference_position] = column + 1;
        }
    }
    if (reference_position != raw_reference_length) {
        return false;
    }

    for (auto& group : groups) {
        group.left_reference = group.minimum_position >
                kCrossAnchorFlankLength
            ? group.minimum_position - kCrossAnchorFlankLength
            : 0;
        const uint64_t requested_right =
            static_cast<uint64_t>(group.maximum_position) +
            kCrossAnchorFlankLength + 1;
        group.right_reference = static_cast<uint32_t>(std::min<uint64_t>(
            raw_reference_length, requested_right));
    }

    std::sort(groups.begin(), groups.end(),
        [](const auto& left, const auto& right) {
            if (left.species_count != right.species_count) {
                return left.species_count > right.species_count;
            }
            if (left.edges.size() != right.edges.size()) {
                return left.edges.size() > right.edges.size();
            }
            if (left.shared_bases != right.shared_bases) {
                return left.shared_bases > right.shared_bases;
            }
            if (left.average_identity != right.average_identity) {
                return left.average_identity > right.average_identity;
            }
            return left.minimum_position < right.minimum_position;
        });

    std::vector<CrossAnchorGroup> selected;
    for (auto& group : groups) {
        const bool overlaps = std::any_of(
            selected.begin(), selected.end(),
            [&](const auto& existing) {
                return rangesOverlap(group, existing);
            });
        if (!overlaps) {
            selected.push_back(std::move(group));
        }
    }
    std::sort(selected.begin(), selected.end(),
        [](const auto& left, const auto& right) {
            return left.left_reference > right.left_reference;
        });

    bool changed = false;
    for (const auto& group : selected) {
        const size_t start_column =
            reference_boundaries[group.left_reference];
        const size_t end_column =
            group.right_reference == raw_reference_length
            ? aligned_reference.size()
            : reference_boundaries[group.right_reference];
        if (start_column >= end_column) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        std::unordered_map<ChrName, std::string> old_local;
        std::set<ChrName> insertion_species;
        for (const size_t event_index : group.event_indices) {
            if (event_index < events.size()) {
                insertion_species.insert(events[event_index].key);
            }
        }
        std::set<ChrName> repaired_keys = insertion_species;
        repaired_keys.insert(reference_name);
        std::map<ChrName, std::string> raw_alleles;
        bool invalid = false;
        for (const auto& [key, sequence] : sequences) {
            if (end_column > sequence.size()) {
                invalid = true;
                break;
            }
            const std::string local = sequence.substr(
                start_column, end_column - start_column);
            old_local.emplace(key, local);
            if (repaired_keys.find(key) == repaired_keys.end()) {
                continue;
            }
            std::string allele = removeAlignmentGaps(local);
            if (allele.size() > configuration.maximum_window_span) {
                invalid = true;
                break;
            }
            raw_alleles.emplace(key, std::move(allele));
        }
        if (invalid) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        std::unordered_map<ChrName, std::string> repaired_local;
        if (!alignLocalAlleles(
                raw_alleles, insertion_species.size(),
                configuration, repaired_local) ||
            repaired_local.size() != raw_alleles.size()) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        std::unordered_map<ChrName, std::string> new_local;
        if (!mergeRepairedAndLegacyProfiles(
                reference_name, insertion_species,
                repaired_local, old_local, new_local)) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        size_t new_length = 0;
        for (const auto& [key, old_row] : old_local) {
            const auto aligned_it = new_local.find(key);
            if (aligned_it == new_local.end() ||
                removeAlignmentGaps(aligned_it->second) !=
                    removeAlignmentGaps(old_row) ||
                (new_length != 0 &&
                 aligned_it->second.size() != new_length)) {
                invalid = true;
                break;
            }
            new_length = aligned_it->second.size();
        }
        if (invalid || new_length == 0) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        bool strict_improvement = false;
        bool degraded = false;
        for (const auto& edge : group.edges) {
            const auto& left_key = events[edge.left].key;
            const auto& right_key = events[edge.right].key;
            const auto old_metrics = alignedPairMetrics(
                old_local.at(left_key), old_local.at(right_key));
            const auto new_metrics = alignedPairMetrics(
                new_local.at(left_key), new_local.at(right_key));
            const double new_identity = new_metrics.paired_bases == 0
                ? 0.0
                : static_cast<double>(new_metrics.matching_bases) /
                    static_cast<double>(new_metrics.paired_bases);
            const uint32_t left_position =
                events[edge.left].reference_position;
            const uint32_t right_position =
                events[edge.right].reference_position;
            const size_t anchor_distance = left_position > right_position
                ? left_position - right_position
                : right_position - left_position;
            // A d-base reference bridge can force up to d homologous
            // insertion residues onto opposite sides of the bridge while
            // preserving sequence order.  Require the full 10-column gain
            // when possible, but do not make an exactly 10 bp insertion
            // impossible to repair merely because its anchors differ.
            const size_t recoverable_gain = edge.paired_bases >
                    anchor_distance
                ? edge.paired_bases - anchor_distance
                : 1;
            const size_t required_gain = std::min<size_t>(
                kCrossAnchorMinimumInsertion, recoverable_gain);
            degraded = degraded ||
                new_metrics.matching_bases < old_metrics.matching_bases;
            strict_improvement = strict_improvement ||
                (new_metrics.paired_bases >=
                     old_metrics.paired_bases +
                         required_gain &&
                 new_identity >= kCrossAnchorMinimumIdentity);
        }
        if (degraded || !strict_improvement) {
            cross_anchor_counters.fallback.fetch_add(
                1, std::memory_order_relaxed);
            continue;
        }

        for (auto& [key, sequence] : sequences) {
            sequence.replace(
                start_column, end_column - start_column,
                new_local.at(key));
        }
        cross_anchor_counters.accepted.fetch_add(
            1, std::memory_order_relaxed);
        changed = true;
    }
    return changed;
}

}  // namespace

void configureCrossAnchorInsertionRepair(
    const std::string& executable,
    uint_t maximum_window_span) {
    {
        std::lock_guard lock(cross_anchor_configuration_mutex);
        cross_anchor_configuration.executable = executable;
        cross_anchor_configuration.maximum_window_span =
            std::max<uint_t>(1, maximum_window_span);
    }
    {
        std::lock_guard lock(cross_anchor_cache_mutex);
        cross_anchor_msa_cache.clear();
    }
    cross_anchor_counters.blocks_scanned.store(0);
    cross_anchor_counters.cigar_candidates.store(0);
    cross_anchor_counters.short_insertions_skipped.store(0);
    cross_anchor_counters.same_anchor_skipped.store(0);
    cross_anchor_counters.distance_skipped.store(0);
    cross_anchor_counters.similarity_pairs.store(0);
    cross_anchor_counters.ksw2_repairs.store(0);
    cross_anchor_counters.minipoa_calls.store(0);
    cross_anchor_counters.cache_hits.store(0);
    cross_anchor_counters.accepted.store(0);
    cross_anchor_counters.fallback.store(0);
    cross_anchor_counters.nanoseconds.store(0);
}

void logCrossAnchorInsertionRepairStats() {
    constexpr double kNanosecondsPerSecond = 1.0e9;
    spdlog::info(
        "[cross-anchor-insertion] blocks_scanned={} cigar_candidates={} "
        "short_insertions_skipped={} same_anchor_skipped={} "
        "distance_skipped={} similar_pairs={} ksw2_repaired={} "
        "minipoa_calls={} cache_hits={} accepted={} fallback={} "
        "wall_seconds={:.3f}",
        cross_anchor_counters.blocks_scanned.load(),
        cross_anchor_counters.cigar_candidates.load(),
        cross_anchor_counters.short_insertions_skipped.load(),
        cross_anchor_counters.same_anchor_skipped.load(),
        cross_anchor_counters.distance_skipped.load(),
        cross_anchor_counters.similarity_pairs.load(),
        cross_anchor_counters.ksw2_repairs.load(),
        cross_anchor_counters.minipoa_calls.load(),
        cross_anchor_counters.cache_hits.load(),
        cross_anchor_counters.accepted.load(),
        cross_anchor_counters.fallback.load(),
        cross_anchor_counters.nanoseconds.load() /
            kNanosecondsPerSecond);
}

uint_t mergeAlignmentByRef(
    ChrName ref_name,
    std::unordered_map<ChrName, std::string>& seqs,
    const std::unordered_map<ChrName, Cigar_t>& cigars) {
    const auto started = std::chrono::steady_clock::now();
    cross_anchor_counters.blocks_scanned.fetch_add(
        1, std::memory_order_relaxed);
    const auto reference_it = seqs.find(ref_name);
    if (reference_it == seqs.end()) {
        throw std::invalid_argument("mergeAlignmentByRef: ref not found");
    }
    const size_t raw_reference_length = reference_it->second.size();
    const auto events = collectInsertionEvents(ref_name, seqs, cigars);
    auto groups = buildCrossAnchorGroups(events, seqs);

    const uint_t legacy_length = mergeAlignmentByRefCore(
        ref_name, seqs, cigars);
    if (!groups.empty()) {
        CrossAnchorRepairConfiguration configuration;
        {
            std::lock_guard lock(cross_anchor_configuration_mutex);
            configuration = cross_anchor_configuration;
        }
        repairCrossAnchorGroups(
            ref_name, raw_reference_length, events,
            std::move(groups), configuration, seqs);
    }
    const auto finished = std::chrono::steady_clock::now();
    cross_anchor_counters.nanoseconds.fetch_add(
        elapsedNanoseconds(started, finished),
        std::memory_order_relaxed);
    const auto repaired_reference = seqs.find(ref_name);
    return repaired_reference == seqs.end()
        ? legacy_length
        : static_cast<uint_t>(repaired_reference->second.size());
}

AlignCount countAlignedBases(const Cigar_t& cigar) {
    AlignCount cnt;
    for (auto op : cigar) {
        uint32_t len;
        char type;
        intToCigar(op, type, len);
        switch (type) {
        case 'M': // match or mismatch
        case '=': // match
        case 'X': // mismatch
            cnt.ref_bases += len;
            cnt.query_bases += len;
            break;
        case 'I': // insertion wrt ref
            cnt.query_bases += len;
            break;
        case 'D': // deletion wrt ref
            cnt.ref_bases += len;
            break;
            // 视情况处理 clip/skip
        case 'S': // soft clip
            cnt.query_bases += len;
            break;
        case 'H': // hard clip
            // 不计入
            break;
        case 'N': // skipped region in ref
            cnt.ref_bases += len;
            break;
        default:
            break;
        }
    }
    return cnt;
}
