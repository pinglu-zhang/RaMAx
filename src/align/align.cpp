#include "align.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
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
    total_length = mergeAlignmentByRef(ref_name, seqs, cigars);
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
uint_t mergeAlignmentByRef(
    ChrName ref_name,
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

