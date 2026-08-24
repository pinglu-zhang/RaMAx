// RaMAx.cpp: 定义应用程序的入口点。
// 主程序：负责解析命令行参数，初始化配置，执行多基因组比对流程

#include "SeqPro.h"
#include "data_process.h"
#include "dependency_preflight.h"
#include "config.hpp"
#include "index.h"
#include "minipoa_locator.h"
#include "output_spec.hpp"
#include "ramax_version.h"
#include "rare_aligner.h"
#include "sequence_utils.h"
#include "external_msa_runner.h"
#include <cstdlib>

#include <array>
#include <cmath>
#include <regex>
#include <unordered_set>

// ------------------------------------------------------------------
// 通用命令行参数结构体（支持 cereal 序列化）
// 用于控制基因组比对 / 搜索 / 组装相关流程的全局参数
// ------------------------------------------------------------------
struct CommonArgs {

    uint32_t schema_version = 6;

    // ========================
    // 输入 / 输出路径相关参数
    // ========================

    std::filesystem::path input_path = "";      // 输入序列文件路径seqfile
    std::filesystem::path output_path = "";     // 最终输出结果路径
    std::vector<std::filesystem::path> output_paths; // CLI/restart multi-output list
    std::vector<RaMAxOutput::OutputSpec> outputs;     // validated, canonical order
    std::filesystem::path work_dir_path = "";   // 工作目录路径（用于中间文件与缓存）


    // ========================
    // 运行控制参数
    // ========================

    bool restart = false;                       // 是否从已有中间结果重新启动任务
    int thread_num = std::thread::hardware_concurrency(); // 使用的线程数（默认使用所有 CPU 核心）

    MultipleGenomeOutputFormat output_format =
        MultipleGenomeOutputFormat::UNKNOWN;    // 多基因组输出格式（HAL、MAF、PAF、GFA）
    std::string paf_mode = "connected";          // PAF pairing: connected or all
    bool paf_mode_explicit = false;              // CLI-only; intentionally not serialized
    std::string gfa_version = "1.1";             // GFA path encoding: 1.0 P-lines or 1.1 W-lines
    bool gfa_version_explicit = false;            // CLI-only; intentionally not serialized
    std::string gfa_profile = "exact";            // exact audit graph or compact-v2-balanced graph
    bool gfa_profile_explicit = false;            // CLI-only; intentionally not serialized
    bool trust_legacy_cache = false;              // persisted until legacy workdir is discarded
    bool pending_config_update = false;           // runtime-only post-input-validation write

    bool enable_repeat_masking = false;          // 是否启用重复序列遮蔽（Repeat Masking）

    // ========================
    // 搜索 / 算法相关参数
    // ========================

    SearchMode search_mode = ACCURATE_SEARCH;   // 搜索模式（如精确搜索 / 快速搜索）
    uint_t accurate_skip_threshold = 10000;     // accurate 模式长 MUM 跳跃阈值（0=关闭）
    bool allow_MEM = false;                     // 是否允许使用 MEM（Maximal Exact Match）
    bool fast_build = true;                     // 是否启用快速索引构建模式
    SeqPro::Length sampling_interval = 32;      // 索引采样间隔（影响速度与内存）
    uint_t sa_sampling_rate = 1;                // 后缀数组文本位置采样率（1=完整 SA）
    uint_t min_span = 65;                       // 锚点或匹配的最小跨度阈值
    double near_distance = 0.01;                // 首轮 d < 此值使用 wfmash
    double far_distance = 0.02;                 // 预留的远缘物种阈值

    // ========================
    // 基因组切片与锚点参数
    // ========================

    uint_t chunk_size = 10000000;               // 基因组切片长度（bp）
    uint_t overlap_size = 0;                    // 相邻切片之间的重叠长度（bp）
    uint_t min_anchor_length = 20;              // 锚点（anchor）的最小长度
    uint_t max_anchor_frequency = 50;           // 单个锚点在基因组中允许的最大出现次数

    // ========================
    // 日志与输出控制
    // ========================

    std::string log_level = "info";              // 日志级别（debug / info / warn / error）
    bool verbose = false;                       // 是否启用详细输出模式
    bool quiet = false;                         // 是否启用静默模式（尽量减少输出）

    // ========================
    // HAL / 多基因组相关参数
    // ========================

    std::string root_name = "";                  // HAL 文件中的根基因组名称
    std::string ref_name = "";                   // 参考基因组名称
    bool merge_exact_contiguous_blocks = true;   // 每轮合并连续 Block
    uint_t merge_query_gap_max = 100;            // query 正间隔上限；0=仅严格连续
    bool one_round = false;                     // 是否只执行一轮处理流程

    bool realign_single_missing_species = true;
    uint_t species_mismatch_realign_max_span = 3000;
    uint_t species_mismatch_zero_gap_max_span = 200;
    bool repair_structural_breaks = true;
    uint_t structural_break_max_span = 1000;
    bool repair_short_blocks = true;

    // ========================
    // cereal 序列化支持
    // ========================

    template<class Archive>
    void serialize(Archive& ar) {
        ar(
            CEREAL_NVP(schema_version),
            CEREAL_NVP(input_path),
            CEREAL_NVP(output_path),
            CEREAL_NVP(output_paths),
            CEREAL_NVP(work_dir_path),
            CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size),
            CEREAL_NVP(min_anchor_length),
            CEREAL_NVP(max_anchor_frequency),
            CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format),
            CEREAL_NVP(paf_mode),
            CEREAL_NVP(gfa_version),
            CEREAL_NVP(gfa_profile),
            CEREAL_NVP(trust_legacy_cache),
            CEREAL_NVP(enable_repeat_masking),
            CEREAL_NVP(search_mode),
            CEREAL_NVP(accurate_skip_threshold),
            CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build),
            CEREAL_NVP(sampling_interval),
            CEREAL_NVP(sa_sampling_rate),
            CEREAL_NVP(min_span),
            CEREAL_NVP(near_distance),
            CEREAL_NVP(far_distance),
            CEREAL_NVP(log_level),
            CEREAL_NVP(verbose),
            CEREAL_NVP(quiet),
            CEREAL_NVP(root_name),
            CEREAL_NVP(ref_name),
            CEREAL_NVP(one_round),
            CEREAL_NVP(merge_exact_contiguous_blocks),
            CEREAL_NVP(merge_query_gap_max),
            CEREAL_NVP(realign_single_missing_species),
            CEREAL_NVP(species_mismatch_realign_max_span),
            CEREAL_NVP(species_mismatch_zero_gap_max_span),
            CEREAL_NVP(repair_structural_breaks),
            CEREAL_NVP(structural_break_max_span),
            CEREAL_NVP(repair_short_blocks)
        );
    }
};

namespace {
constexpr uint32_t CONFIG_SCHEMA_VERSION = 6;
constexpr uint32_t OUTPUTS_SCHEMA_VERSION = 1;
constexpr const char* OUTPUTS_FILE = "outputs.json";
constexpr const char* INPUT_MANIFEST_FILE = "input_manifest.json";

struct OutputManifest {
    uint32_t schema_version = OUTPUTS_SCHEMA_VERSION;
    std::vector<FilePath> output_paths;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CEREAL_NVP(schema_version), CEREAL_NVP(output_paths));
    }
};

// Exact reader for the configuration written by RaMAx 1.0.6.  Several
// parameters were not serialized in that schema and intentionally retain the
// 1.0.6 restart defaults when migrated.
struct LegacyCommonArgsV1 {
    uint32_t schema_version{1};
    FilePath input_path;
    FilePath output_path;
    FilePath work_dir_path;
    uint_t chunk_size{10000000};
    uint_t overlap_size{0};
    bool restart{false};
    int thread_num{static_cast<int>(std::thread::hardware_concurrency())};
    MultipleGenomeOutputFormat output_format{MultipleGenomeOutputFormat::UNKNOWN};
    bool enable_repeat_masking{false};
    SearchMode search_mode{ACCURATE_SEARCH};
    bool allow_MEM{false};
    bool fast_build{true};
    SeqPro::Length sampling_interval{32};
    uint_t min_span{65};
    std::string log_level{"info"};
    bool verbose{false};
    bool quiet{false};
    std::string root_name;
    std::string ref_name;
    bool one_round{false};
    bool merge_exact_contiguous_blocks{true};
    uint_t merge_query_gap_max{100};
    bool realign_single_missing_species{true};
    uint_t species_mismatch_realign_max_span{3000};
    uint_t species_mismatch_zero_gap_max_span{200};
    bool repair_structural_breaks{true};
    uint_t structural_break_max_span{1000};
    bool repair_short_blocks{true};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            CEREAL_NVP(schema_version), CEREAL_NVP(input_path),
            CEREAL_NVP(output_path), CEREAL_NVP(work_dir_path),
            CEREAL_NVP(chunk_size), CEREAL_NVP(overlap_size),
            CEREAL_NVP(restart), CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format), CEREAL_NVP(enable_repeat_masking),
            CEREAL_NVP(search_mode), CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build), CEREAL_NVP(sampling_interval),
            CEREAL_NVP(min_span), CEREAL_NVP(log_level),
            CEREAL_NVP(verbose), CEREAL_NVP(quiet),
            CEREAL_NVP(root_name), CEREAL_NVP(ref_name),
            CEREAL_NVP(one_round),
            CEREAL_NVP(merge_exact_contiguous_blocks),
            CEREAL_NVP(merge_query_gap_max),
            CEREAL_NVP(realign_single_missing_species),
            CEREAL_NVP(species_mismatch_realign_max_span),
            CEREAL_NVP(species_mismatch_zero_gap_max_span),
            CEREAL_NVP(repair_structural_breaks),
            CEREAL_NVP(structural_break_max_span),
            CEREAL_NVP(repair_short_blocks));
    }
};

// Exact reader for schema 2, before Mash-distance routing thresholds became
// persistent restart parameters.
struct LegacyCommonArgsV2 {
    uint32_t schema_version{2};
    FilePath input_path;
    FilePath output_path;
    std::vector<FilePath> output_paths;
    FilePath work_dir_path;
    uint_t chunk_size{10000000};
    uint_t overlap_size{0};
    uint_t min_anchor_length{20};
    uint_t max_anchor_frequency{50};
    int thread_num{static_cast<int>(std::thread::hardware_concurrency())};
    MultipleGenomeOutputFormat output_format{MultipleGenomeOutputFormat::UNKNOWN};
    std::string paf_mode{"connected"};
    bool trust_legacy_cache{false};
    bool enable_repeat_masking{false};
    SearchMode search_mode{ACCURATE_SEARCH};
    uint_t accurate_skip_threshold{10000};
    bool allow_MEM{false};
    bool fast_build{true};
    SeqPro::Length sampling_interval{32};
    uint_t min_span{65};
    std::string log_level{"info"};
    bool verbose{false};
    bool quiet{false};
    std::string root_name;
    std::string ref_name;
    bool one_round{false};
    bool merge_exact_contiguous_blocks{true};
    uint_t merge_query_gap_max{100};
    bool realign_single_missing_species{true};
    uint_t species_mismatch_realign_max_span{3000};
    uint_t species_mismatch_zero_gap_max_span{200};
    bool repair_structural_breaks{true};
    uint_t structural_break_max_span{1000};
    bool repair_short_blocks{true};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            CEREAL_NVP(schema_version), CEREAL_NVP(input_path),
            CEREAL_NVP(output_path), CEREAL_NVP(output_paths),
            CEREAL_NVP(work_dir_path), CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size), CEREAL_NVP(min_anchor_length),
            CEREAL_NVP(max_anchor_frequency), CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format), CEREAL_NVP(paf_mode),
            CEREAL_NVP(trust_legacy_cache),
            CEREAL_NVP(enable_repeat_masking), CEREAL_NVP(search_mode),
            CEREAL_NVP(accurate_skip_threshold), CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build), CEREAL_NVP(sampling_interval),
            CEREAL_NVP(min_span), CEREAL_NVP(log_level),
            CEREAL_NVP(verbose), CEREAL_NVP(quiet),
            CEREAL_NVP(root_name), CEREAL_NVP(ref_name),
            CEREAL_NVP(one_round),
            CEREAL_NVP(merge_exact_contiguous_blocks),
            CEREAL_NVP(merge_query_gap_max),
            CEREAL_NVP(realign_single_missing_species),
            CEREAL_NVP(species_mismatch_realign_max_span),
            CEREAL_NVP(species_mismatch_zero_gap_max_span),
            CEREAL_NVP(repair_structural_breaks),
            CEREAL_NVP(structural_break_max_span),
            CEREAL_NVP(repair_short_blocks));
    }
};

// Exact reader for schema 3, before the selectable GFA path encoding became
// a persistent restart parameter.
struct LegacyCommonArgsV3 {
    uint32_t schema_version{3};
    FilePath input_path;
    FilePath output_path;
    std::vector<FilePath> output_paths;
    FilePath work_dir_path;
    uint_t chunk_size{10000000};
    uint_t overlap_size{0};
    uint_t min_anchor_length{20};
    uint_t max_anchor_frequency{50};
    int thread_num{static_cast<int>(std::thread::hardware_concurrency())};
    MultipleGenomeOutputFormat output_format{MultipleGenomeOutputFormat::UNKNOWN};
    std::string paf_mode{"connected"};
    bool trust_legacy_cache{false};
    bool enable_repeat_masking{false};
    SearchMode search_mode{ACCURATE_SEARCH};
    uint_t accurate_skip_threshold{10000};
    bool allow_MEM{false};
    bool fast_build{true};
    SeqPro::Length sampling_interval{32};
    uint_t min_span{65};
    double near_distance{0.01};
    double far_distance{0.02};
    std::string log_level{"info"};
    bool verbose{false};
    bool quiet{false};
    std::string root_name;
    std::string ref_name;
    bool one_round{false};
    bool merge_exact_contiguous_blocks{true};
    uint_t merge_query_gap_max{100};
    bool realign_single_missing_species{true};
    uint_t species_mismatch_realign_max_span{3000};
    uint_t species_mismatch_zero_gap_max_span{200};
    bool repair_structural_breaks{true};
    uint_t structural_break_max_span{1000};
    bool repair_short_blocks{true};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            CEREAL_NVP(schema_version), CEREAL_NVP(input_path),
            CEREAL_NVP(output_path), CEREAL_NVP(output_paths),
            CEREAL_NVP(work_dir_path), CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size), CEREAL_NVP(min_anchor_length),
            CEREAL_NVP(max_anchor_frequency), CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format), CEREAL_NVP(paf_mode),
            CEREAL_NVP(trust_legacy_cache),
            CEREAL_NVP(enable_repeat_masking), CEREAL_NVP(search_mode),
            CEREAL_NVP(accurate_skip_threshold), CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build), CEREAL_NVP(sampling_interval),
            CEREAL_NVP(min_span), CEREAL_NVP(near_distance),
            CEREAL_NVP(far_distance), CEREAL_NVP(log_level),
            CEREAL_NVP(verbose), CEREAL_NVP(quiet),
            CEREAL_NVP(root_name), CEREAL_NVP(ref_name),
            CEREAL_NVP(one_round),
            CEREAL_NVP(merge_exact_contiguous_blocks),
            CEREAL_NVP(merge_query_gap_max),
            CEREAL_NVP(realign_single_missing_species),
            CEREAL_NVP(species_mismatch_realign_max_span),
            CEREAL_NVP(species_mismatch_zero_gap_max_span),
            CEREAL_NVP(repair_structural_breaks),
            CEREAL_NVP(structural_break_max_span),
            CEREAL_NVP(repair_short_blocks));
    }
};

// Exact reader for schema 4, before the selectable GFA graph profile became
// a persistent restart parameter.
struct LegacyCommonArgsV4 {
    uint32_t schema_version{4};
    FilePath input_path;
    FilePath output_path;
    std::vector<FilePath> output_paths;
    FilePath work_dir_path;
    uint_t chunk_size{10000000};
    uint_t overlap_size{0};
    uint_t min_anchor_length{20};
    uint_t max_anchor_frequency{50};
    int thread_num{static_cast<int>(std::thread::hardware_concurrency())};
    MultipleGenomeOutputFormat output_format{MultipleGenomeOutputFormat::UNKNOWN};
    std::string paf_mode{"connected"};
    std::string gfa_version{"1.1"};
    bool trust_legacy_cache{false};
    bool enable_repeat_masking{false};
    SearchMode search_mode{ACCURATE_SEARCH};
    uint_t accurate_skip_threshold{10000};
    bool allow_MEM{false};
    bool fast_build{true};
    SeqPro::Length sampling_interval{32};
    uint_t min_span{65};
    double near_distance{0.01};
    double far_distance{0.02};
    std::string log_level{"info"};
    bool verbose{false};
    bool quiet{false};
    std::string root_name;
    std::string ref_name;
    bool one_round{false};
    bool merge_exact_contiguous_blocks{true};
    uint_t merge_query_gap_max{100};
    bool realign_single_missing_species{true};
    uint_t species_mismatch_realign_max_span{3000};
    uint_t species_mismatch_zero_gap_max_span{200};
    bool repair_structural_breaks{true};
    uint_t structural_break_max_span{1000};
    bool repair_short_blocks{true};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            CEREAL_NVP(schema_version), CEREAL_NVP(input_path),
            CEREAL_NVP(output_path), CEREAL_NVP(output_paths),
            CEREAL_NVP(work_dir_path), CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size), CEREAL_NVP(min_anchor_length),
            CEREAL_NVP(max_anchor_frequency), CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format), CEREAL_NVP(paf_mode),
            CEREAL_NVP(gfa_version), CEREAL_NVP(trust_legacy_cache),
            CEREAL_NVP(enable_repeat_masking), CEREAL_NVP(search_mode),
            CEREAL_NVP(accurate_skip_threshold), CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build), CEREAL_NVP(sampling_interval),
            CEREAL_NVP(min_span), CEREAL_NVP(near_distance),
            CEREAL_NVP(far_distance), CEREAL_NVP(log_level),
            CEREAL_NVP(verbose), CEREAL_NVP(quiet),
            CEREAL_NVP(root_name), CEREAL_NVP(ref_name),
            CEREAL_NVP(one_round),
            CEREAL_NVP(merge_exact_contiguous_blocks),
            CEREAL_NVP(merge_query_gap_max),
            CEREAL_NVP(realign_single_missing_species),
            CEREAL_NVP(species_mismatch_realign_max_span),
            CEREAL_NVP(species_mismatch_zero_gap_max_span),
            CEREAL_NVP(repair_structural_breaks),
            CEREAL_NVP(structural_break_max_span),
            CEREAL_NVP(repair_short_blocks));
    }
};


// Exact reader for schema 5, before suffix-array text-position sampling
// became a persistent restart parameter.
struct LegacyCommonArgsV5 {
    uint32_t schema_version{5};
    FilePath input_path;
    FilePath output_path;
    std::vector<FilePath> output_paths;
    FilePath work_dir_path;
    uint_t chunk_size{10000000};
    uint_t overlap_size{0};
    uint_t min_anchor_length{20};
    uint_t max_anchor_frequency{50};
    int thread_num{static_cast<int>(std::thread::hardware_concurrency())};
    MultipleGenomeOutputFormat output_format{MultipleGenomeOutputFormat::UNKNOWN};
    std::string paf_mode{"connected"};
    std::string gfa_version{"1.1"};
    std::string gfa_profile{"exact"};
    bool trust_legacy_cache{false};
    bool enable_repeat_masking{false};
    SearchMode search_mode{ACCURATE_SEARCH};
    uint_t accurate_skip_threshold{10000};
    bool allow_MEM{false};
    bool fast_build{true};
    SeqPro::Length sampling_interval{32};
    uint_t min_span{65};
    double near_distance{0.01};
    double far_distance{0.02};
    std::string log_level{"info"};
    bool verbose{false};
    bool quiet{false};
    std::string root_name;
    std::string ref_name;
    bool one_round{false};
    bool merge_exact_contiguous_blocks{true};
    uint_t merge_query_gap_max{100};
    bool realign_single_missing_species{true};
    uint_t species_mismatch_realign_max_span{3000};
    uint_t species_mismatch_zero_gap_max_span{200};
    bool repair_structural_breaks{true};
    uint_t structural_break_max_span{1000};
    bool repair_short_blocks{true};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            CEREAL_NVP(schema_version), CEREAL_NVP(input_path),
            CEREAL_NVP(output_path), CEREAL_NVP(output_paths),
            CEREAL_NVP(work_dir_path), CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size), CEREAL_NVP(min_anchor_length),
            CEREAL_NVP(max_anchor_frequency), CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format), CEREAL_NVP(paf_mode),
            CEREAL_NVP(gfa_version), CEREAL_NVP(gfa_profile),
            CEREAL_NVP(trust_legacy_cache),
            CEREAL_NVP(enable_repeat_masking), CEREAL_NVP(search_mode),
            CEREAL_NVP(accurate_skip_threshold), CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build), CEREAL_NVP(sampling_interval),
            CEREAL_NVP(min_span), CEREAL_NVP(near_distance),
            CEREAL_NVP(far_distance), CEREAL_NVP(log_level),
            CEREAL_NVP(verbose), CEREAL_NVP(quiet),
            CEREAL_NVP(root_name), CEREAL_NVP(ref_name),
            CEREAL_NVP(one_round),
            CEREAL_NVP(merge_exact_contiguous_blocks),
            CEREAL_NVP(merge_query_gap_max),
            CEREAL_NVP(realign_single_missing_species),
            CEREAL_NVP(species_mismatch_realign_max_span),
            CEREAL_NVP(species_mismatch_zero_gap_max_span),
            CEREAL_NVP(repair_structural_breaks),
            CEREAL_NVP(structural_break_max_span),
            CEREAL_NVP(repair_short_blocks));
    }
};

struct InputIdentityRecord {
    std::string species;
    std::string source;
    bool source_is_url{false};
    uint64_t source_size{0};
    int64_t source_mtime{0};

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CEREAL_NVP(species), CEREAL_NVP(source),
                CEREAL_NVP(source_is_url), CEREAL_NVP(source_size),
                CEREAL_NVP(source_mtime));
    }

    bool operator==(const InputIdentityRecord&) const = default;
};

struct InputManifest {
    uint32_t schema_version{1};
    std::string seqfile;
    uint64_t seqfile_size{0};
    int64_t seqfile_mtime{0};
    std::vector<InputIdentityRecord> inputs;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(CEREAL_NVP(schema_version), CEREAL_NVP(seqfile),
                CEREAL_NVP(seqfile_size), CEREAL_NVP(seqfile_mtime),
                CEREAL_NVP(inputs));
    }
};

struct RestartOverrides {
    CommonArgs values;
    std::unordered_set<std::string> specified;

    bool has(const std::string& name) const {
        return specified.contains(name);
    }
};

uint32_t detectConfigSchema(const FilePath& config_path) {
    std::ifstream input(config_path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Failed to open restart configuration: " +
                                 config_path.string());
    }
    const std::string text((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());
    const std::regex pattern("\\\"schema_version\\\"\\s*:\\s*([0-9]+)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        throw std::runtime_error(
            "Restart configuration has no schema_version: " +
            config_path.string());
    }
    return static_cast<uint32_t>(std::stoul(match[1].str()));
}

void saveEffectiveConfig(const CommonArgs& args) {
    FilePath config_path = args.work_dir_path / CONFIG_FILE;
    FilePath partial = config_path;
    partial += ".partial";
    RaMAxCache::removeIfPresent(partial);
    CommonArgs saved = args;
    saved.schema_version = CONFIG_SCHEMA_VERSION;
    saved.restart = false;
    saved.pending_config_update = false;
    try {
        std::ofstream output(partial, std::ios::binary | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Failed to save configuration: " +
                                     partial.string());
        }
        {
            cereal::JSONOutputArchive archive(output);
            archive(cereal::make_nvp("common_args", saved));
        }
        output.close();
        RaMAxCache::publishFile(partial, config_path);
    } catch (...) {
        RaMAxCache::removeIfPresent(partial);
        throw;
    }
}

CommonArgs convertLegacyConfig(const LegacyCommonArgsV1& legacy) {
    CommonArgs args;
    args.schema_version = CONFIG_SCHEMA_VERSION;
    args.input_path = legacy.input_path;
    args.output_path = legacy.output_path;
    args.work_dir_path = legacy.work_dir_path;
    args.chunk_size = legacy.chunk_size;
    args.overlap_size = legacy.overlap_size;
    args.thread_num = legacy.thread_num;
    args.output_format = legacy.output_format;
    args.enable_repeat_masking = legacy.enable_repeat_masking;
    args.search_mode = legacy.search_mode;
    args.allow_MEM = legacy.allow_MEM;
    args.fast_build = legacy.fast_build;
    args.sampling_interval = legacy.sampling_interval;
    args.min_span = legacy.min_span;
    args.log_level = legacy.log_level;
    args.verbose = legacy.verbose;
    args.quiet = legacy.quiet;
    args.root_name = legacy.root_name;
    args.ref_name = legacy.ref_name;
    args.one_round = legacy.one_round;
    args.merge_exact_contiguous_blocks = legacy.merge_exact_contiguous_blocks;
    args.merge_query_gap_max = legacy.merge_query_gap_max;
    args.realign_single_missing_species = legacy.realign_single_missing_species;
    args.species_mismatch_realign_max_span =
        legacy.species_mismatch_realign_max_span;
    args.species_mismatch_zero_gap_max_span =
        legacy.species_mismatch_zero_gap_max_span;
    args.repair_structural_breaks = legacy.repair_structural_breaks;
    args.structural_break_max_span = legacy.structural_break_max_span;
    args.repair_short_blocks = legacy.repair_short_blocks;
    // Fields absent from schema 1 retain the exact 1.0.6 restart defaults.
    args.min_anchor_length = 20;
    args.max_anchor_frequency = 50;
    args.accurate_skip_threshold = 10000;
    args.paf_mode = "connected";
    args.trust_legacy_cache = true;
    return args;
}

CommonArgs convertSchema2Config(const LegacyCommonArgsV2& legacy) {
    CommonArgs args;
    args.schema_version = CONFIG_SCHEMA_VERSION;
    args.input_path = legacy.input_path;
    args.output_path = legacy.output_path;
    args.output_paths = legacy.output_paths;
    args.work_dir_path = legacy.work_dir_path;
    args.chunk_size = legacy.chunk_size;
    args.overlap_size = legacy.overlap_size;
    args.min_anchor_length = legacy.min_anchor_length;
    args.max_anchor_frequency = legacy.max_anchor_frequency;
    args.thread_num = legacy.thread_num;
    args.output_format = legacy.output_format;
    args.paf_mode = legacy.paf_mode;
    args.trust_legacy_cache = legacy.trust_legacy_cache;
    args.enable_repeat_masking = legacy.enable_repeat_masking;
    args.search_mode = legacy.search_mode;
    args.accurate_skip_threshold = legacy.accurate_skip_threshold;
    args.allow_MEM = legacy.allow_MEM;
    args.fast_build = legacy.fast_build;
    args.sampling_interval = legacy.sampling_interval;
    args.min_span = legacy.min_span;
    args.log_level = legacy.log_level;
    args.verbose = legacy.verbose;
    args.quiet = legacy.quiet;
    args.root_name = legacy.root_name;
    args.ref_name = legacy.ref_name;
    args.one_round = legacy.one_round;
    args.merge_exact_contiguous_blocks = legacy.merge_exact_contiguous_blocks;
    args.merge_query_gap_max = legacy.merge_query_gap_max;
    args.realign_single_missing_species = legacy.realign_single_missing_species;
    args.species_mismatch_realign_max_span =
        legacy.species_mismatch_realign_max_span;
    args.species_mismatch_zero_gap_max_span =
        legacy.species_mismatch_zero_gap_max_span;
    args.repair_structural_breaks = legacy.repair_structural_breaks;
    args.structural_break_max_span = legacy.structural_break_max_span;
    args.repair_short_blocks = legacy.repair_short_blocks;
    args.near_distance = 0.01;
    args.far_distance = 0.02;
    return args;
}

CommonArgs convertSchema3Config(const LegacyCommonArgsV3& legacy) {
    CommonArgs args;
    args.schema_version = CONFIG_SCHEMA_VERSION;
    args.input_path = legacy.input_path;
    args.output_path = legacy.output_path;
    args.output_paths = legacy.output_paths;
    args.work_dir_path = legacy.work_dir_path;
    args.chunk_size = legacy.chunk_size;
    args.overlap_size = legacy.overlap_size;
    args.min_anchor_length = legacy.min_anchor_length;
    args.max_anchor_frequency = legacy.max_anchor_frequency;
    args.thread_num = legacy.thread_num;
    args.output_format = legacy.output_format;
    args.paf_mode = legacy.paf_mode;
    args.gfa_version = "1.1";
    args.trust_legacy_cache = legacy.trust_legacy_cache;
    args.enable_repeat_masking = legacy.enable_repeat_masking;
    args.search_mode = legacy.search_mode;
    args.accurate_skip_threshold = legacy.accurate_skip_threshold;
    args.allow_MEM = legacy.allow_MEM;
    args.fast_build = legacy.fast_build;
    args.sampling_interval = legacy.sampling_interval;
    args.min_span = legacy.min_span;
    args.near_distance = legacy.near_distance;
    args.far_distance = legacy.far_distance;
    args.log_level = legacy.log_level;
    args.verbose = legacy.verbose;
    args.quiet = legacy.quiet;
    args.root_name = legacy.root_name;
    args.ref_name = legacy.ref_name;
    args.one_round = legacy.one_round;
    args.merge_exact_contiguous_blocks = legacy.merge_exact_contiguous_blocks;
    args.merge_query_gap_max = legacy.merge_query_gap_max;
    args.realign_single_missing_species = legacy.realign_single_missing_species;
    args.species_mismatch_realign_max_span =
        legacy.species_mismatch_realign_max_span;
    args.species_mismatch_zero_gap_max_span =
        legacy.species_mismatch_zero_gap_max_span;
    args.repair_structural_breaks = legacy.repair_structural_breaks;
    args.structural_break_max_span = legacy.structural_break_max_span;
    args.repair_short_blocks = legacy.repair_short_blocks;
    return args;
}

CommonArgs convertSchema4Config(const LegacyCommonArgsV4& legacy) {
    CommonArgs args;
    args.schema_version = CONFIG_SCHEMA_VERSION;
    args.input_path = legacy.input_path;
    args.output_path = legacy.output_path;
    args.output_paths = legacy.output_paths;
    args.work_dir_path = legacy.work_dir_path;
    args.chunk_size = legacy.chunk_size;
    args.overlap_size = legacy.overlap_size;
    args.min_anchor_length = legacy.min_anchor_length;
    args.max_anchor_frequency = legacy.max_anchor_frequency;
    args.thread_num = legacy.thread_num;
    args.output_format = legacy.output_format;
    args.paf_mode = legacy.paf_mode;
    args.gfa_version = legacy.gfa_version;
    args.gfa_profile = "exact";
    args.trust_legacy_cache = legacy.trust_legacy_cache;
    args.enable_repeat_masking = legacy.enable_repeat_masking;
    args.search_mode = legacy.search_mode;
    args.accurate_skip_threshold = legacy.accurate_skip_threshold;
    args.allow_MEM = legacy.allow_MEM;
    args.fast_build = legacy.fast_build;
    args.sampling_interval = legacy.sampling_interval;
    args.min_span = legacy.min_span;
    args.near_distance = legacy.near_distance;
    args.far_distance = legacy.far_distance;
    args.log_level = legacy.log_level;
    args.verbose = legacy.verbose;
    args.quiet = legacy.quiet;
    args.root_name = legacy.root_name;
    args.ref_name = legacy.ref_name;
    args.one_round = legacy.one_round;
    args.merge_exact_contiguous_blocks = legacy.merge_exact_contiguous_blocks;
    args.merge_query_gap_max = legacy.merge_query_gap_max;
    args.realign_single_missing_species = legacy.realign_single_missing_species;
    args.species_mismatch_realign_max_span =
        legacy.species_mismatch_realign_max_span;
    args.species_mismatch_zero_gap_max_span =
        legacy.species_mismatch_zero_gap_max_span;
    args.repair_structural_breaks = legacy.repair_structural_breaks;
    args.structural_break_max_span = legacy.structural_break_max_span;
    args.repair_short_blocks = legacy.repair_short_blocks;
    return args;
}


CommonArgs convertSchema5Config(const LegacyCommonArgsV5& legacy) {
    CommonArgs args;
    args.schema_version = CONFIG_SCHEMA_VERSION;
    args.input_path = legacy.input_path;
    args.output_path = legacy.output_path;
    args.output_paths = legacy.output_paths;
    args.work_dir_path = legacy.work_dir_path;
    args.chunk_size = legacy.chunk_size;
    args.overlap_size = legacy.overlap_size;
    args.min_anchor_length = legacy.min_anchor_length;
    args.max_anchor_frequency = legacy.max_anchor_frequency;
    args.thread_num = legacy.thread_num;
    args.output_format = legacy.output_format;
    args.paf_mode = legacy.paf_mode;
    args.gfa_version = legacy.gfa_version;
    args.gfa_profile = legacy.gfa_profile;
    args.trust_legacy_cache = legacy.trust_legacy_cache;
    args.enable_repeat_masking = legacy.enable_repeat_masking;
    args.search_mode = legacy.search_mode;
    args.accurate_skip_threshold = legacy.accurate_skip_threshold;
    args.allow_MEM = legacy.allow_MEM;
    args.fast_build = legacy.fast_build;
    args.sampling_interval = legacy.sampling_interval;
    args.sa_sampling_rate = 1;
    args.min_span = legacy.min_span;
    args.near_distance = legacy.near_distance;
    args.far_distance = legacy.far_distance;
    args.log_level = legacy.log_level;
    args.verbose = legacy.verbose;
    args.quiet = legacy.quiet;
    args.root_name = legacy.root_name;
    args.ref_name = legacy.ref_name;
    args.one_round = legacy.one_round;
    args.merge_exact_contiguous_blocks = legacy.merge_exact_contiguous_blocks;
    args.merge_query_gap_max = legacy.merge_query_gap_max;
    args.realign_single_missing_species = legacy.realign_single_missing_species;
    args.species_mismatch_realign_max_span =
        legacy.species_mismatch_realign_max_span;
    args.species_mismatch_zero_gap_max_span =
        legacy.species_mismatch_zero_gap_max_span;
    args.repair_structural_breaks = legacy.repair_structural_breaks;
    args.structural_break_max_span = legacy.structural_break_max_span;
    args.repair_short_blocks = legacy.repair_short_blocks;
    return args;
}

void validateGfaVersion(const CommonArgs& args) {
    if (args.gfa_version != "1.0" && args.gfa_version != "1.1") {
        throw std::runtime_error(
            "--gfa-version must be exactly 1.0 or 1.1");
    }
}

void validateGfaProfile(const CommonArgs& args) {
    if (args.gfa_profile != "exact" && args.gfa_profile != "compact") {
        throw std::runtime_error(
            "--gfa-profile must be exactly exact or compact");
    }
}

void validateDistanceThresholds(const CommonArgs& args) {
    if (!std::isfinite(args.near_distance) ||
        !std::isfinite(args.far_distance) ||
        args.near_distance < 0.0 ||
        args.near_distance >= args.far_distance ||
        args.far_distance > 1.0) {
        throw std::runtime_error(
            "Distance thresholds must satisfy 0 <= --near-distance < "
            "--far-distance <= 1");
    }
}


void validateSaSamplingRate(const CommonArgs& args) {
    if (args.sa_sampling_rate != 1) {
        throw std::runtime_error(
            "--sa-sampling-rate must be exactly 1; the suffix-array backend "
            "stores complete SA/ISA/LCP arrays");
    }
}

RestartOverrides captureRestartOverrides(const CLI::App& app,
                                         const CommonArgs& values) {
    RestartOverrides overrides;
    overrides.values = values;
    constexpr std::array<const char*, 32> names{
        "--output", "--paf-mode", "--gfa-version", "--gfa-profile", "--chunk_size", "--root", "--ref",
        "--overlap_size", "--min_anchor_length", "--max_anchor_frequency",
        "--search-mode", "--accurate-skip-threshold", "--allow-mem",
        "--one-round", "--optimize-blocks", "--merge-blocks", "--merge-gap",
        "--realign-missing", "--realign-span", "--zero-gap-span",
        "--repair-breaks", "--break-span", "--merge-short-blocks",
        "--slow-build", "--sampling-interval", "--sa-sampling-rate",
        "--min-span", "--threads",
        "--log-level", "--verbose", "--near-distance", "--far-distance"
    };
    for (const char* name : names) {
        if (app.count(name) != 0) overrides.specified.emplace(name);
    }
    if (app.count("--quiet") != 0) overrides.specified.emplace("--quiet");
    return overrides;
}

void applyRestartOverrides(CommonArgs& args,
                           const RestartOverrides& overrides) {
    const CommonArgs& value = overrides.values;
    if (overrides.has("--output")) args.output_paths = value.output_paths;
    if (overrides.has("--paf-mode")) args.paf_mode = value.paf_mode;
    if (overrides.has("--gfa-version")) args.gfa_version = value.gfa_version;
    if (overrides.has("--gfa-profile")) args.gfa_profile = value.gfa_profile;
    if (overrides.has("--chunk_size")) args.chunk_size = value.chunk_size;
    if (overrides.has("--root")) args.root_name = value.root_name;
    if (overrides.has("--ref")) args.ref_name = value.ref_name;
    if (overrides.has("--overlap_size")) args.overlap_size = value.overlap_size;
    if (overrides.has("--min_anchor_length"))
        args.min_anchor_length = value.min_anchor_length;
    if (overrides.has("--max_anchor_frequency"))
        args.max_anchor_frequency = value.max_anchor_frequency;
    if (overrides.has("--search-mode")) args.search_mode = value.search_mode;
    if (overrides.has("--accurate-skip-threshold"))
        args.accurate_skip_threshold = value.accurate_skip_threshold;
    if (overrides.has("--allow-mem")) args.allow_MEM = true;
    if (overrides.has("--one-round")) args.one_round = true;
    if (overrides.has("--merge-blocks"))
        args.merge_exact_contiguous_blocks = true;
    if (overrides.has("--merge-gap"))
        args.merge_query_gap_max = value.merge_query_gap_max;
    if (overrides.has("--realign-missing"))
        args.realign_single_missing_species = true;
    if (overrides.has("--realign-span"))
        args.species_mismatch_realign_max_span =
            value.species_mismatch_realign_max_span;
    if (overrides.has("--zero-gap-span"))
        args.species_mismatch_zero_gap_max_span =
            value.species_mismatch_zero_gap_max_span;
    if (overrides.has("--repair-breaks")) args.repair_structural_breaks = true;
    if (overrides.has("--break-span"))
        args.structural_break_max_span = value.structural_break_max_span;
    if (overrides.has("--merge-short-blocks")) args.repair_short_blocks = true;
    if (overrides.has("--slow-build")) args.fast_build = false;
    if (overrides.has("--sampling-interval"))
        args.sampling_interval = value.sampling_interval;
    if (overrides.has("--sa-sampling-rate"))
        args.sa_sampling_rate = value.sa_sampling_rate;
    if (overrides.has("--min-span")) args.min_span = value.min_span;
    if (overrides.has("--near-distance"))
        args.near_distance = value.near_distance;
    if (overrides.has("--far-distance"))
        args.far_distance = value.far_distance;
    if (overrides.has("--threads")) args.thread_num = value.thread_num;
    if (overrides.has("--log-level")) args.log_level = value.log_level;
    if (overrides.has("--verbose")) {
        args.verbose = true;
        args.quiet = false;
    }
    if (overrides.has("--quiet")) {
        args.quiet = true;
        args.verbose = false;
    }
    args.paf_mode_explicit = overrides.has("--paf-mode");
    args.gfa_version_explicit = overrides.has("--gfa-version");
    args.gfa_profile_explicit = overrides.has("--gfa-profile");
}

InputManifest makeInputManifest(const CommonArgs& args,
                                const SpeciesPathMap& all_species) {
    InputManifest manifest;
    manifest.seqfile = args.input_path.string();
    const auto seqfile_metadata = RaMAxCache::fileMetadata(args.input_path);
    manifest.seqfile_size = seqfile_metadata.size;
    manifest.seqfile_mtime = seqfile_metadata.mtime;
    manifest.inputs.reserve(all_species.size());
    for (const auto& [species, path] : all_species) {
        InputIdentityRecord record;
        record.species = species;
        record.source = path.string();
        record.source_is_url = isUrl(path.string());
        if (!record.source_is_url) {
            const auto metadata = RaMAxCache::fileMetadata(path);
            record.source_size = metadata.size;
            record.source_mtime = metadata.mtime;
        }
        manifest.inputs.push_back(std::move(record));
    }
    std::sort(manifest.inputs.begin(), manifest.inputs.end(),
              [](const auto& left, const auto& right) {
                  return left.species < right.species;
              });
    return manifest;
}

void saveInputManifest(const FilePath& path, const InputManifest& manifest) {
    FilePath partial = path;
    partial += ".partial";
    RaMAxCache::removeIfPresent(partial);
    try {
        std::ofstream output(partial, std::ios::binary | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("Failed to save input manifest: " +
                                     partial.string());
        }
        {
            cereal::JSONOutputArchive archive(output);
            archive(cereal::make_nvp("input_manifest", manifest));
        }
        output.close();
        RaMAxCache::publishFile(partial, path);
    } catch (...) {
        RaMAxCache::removeIfPresent(partial);
        throw;
    }
}

void validateOrCreateInputManifest(const CommonArgs& args,
                                   const SpeciesPathMap& all_species) {
    const FilePath path = args.work_dir_path / INPUT_MANIFEST_FILE;
    const InputManifest current = makeInputManifest(args, all_species);
    if (!std::filesystem::exists(path)) {
        if (args.restart) {
            spdlog::warn(
                "Input manifest is absent; trusting current schema-1 input identity and migrating it");
        }
        saveInputManifest(path, current);
        return;
    }

    InputManifest saved;
    try {
        std::ifstream input(path);
        cereal::JSONInputArchive archive(input);
        archive(cereal::make_nvp("input_manifest", saved));
    } catch (const std::exception& error) {
        throw std::runtime_error("Invalid input manifest " + path.string() +
                                 ": " + error.what());
    }
    if (saved.schema_version != current.schema_version ||
        saved.seqfile != current.seqfile ||
        saved.seqfile_size != current.seqfile_size ||
        saved.seqfile_mtime != current.seqfile_mtime ||
        saved.inputs != current.inputs) {
        throw std::runtime_error(
            "Restart input identity changed; use a new work directory instead of reusing " +
            args.work_dir_path.string());
    }
}

bool hasReusableRawSnapshot(const CommonArgs& args,
                            const SpeciesName& species,
                            const FilePath& source) {
    const FilePath raw = args.work_dir_path / DATA_DIR / RAW_DATA_DIR /
        (species + getFileExtension(source));
    const FilePath marker = RaMAxCache::completionMarkerPath(raw);
    if (RaMAxCache::markerMatches(
            marker, "raw-fasta", 1, source.string(), true, source, raw)) {
        return true;
    }
    return args.trust_legacy_cache &&
           std::filesystem::is_regular_file(raw) &&
           !std::filesystem::exists(marker);
}

void archivePreviousLog(const FilePath& work_dir) {
    const FilePath current = work_dir / LOGGER_FILE;
    if (!std::filesystem::is_regular_file(current)) return;
    for (uint32_t index = 1; ; ++index) {
        const FilePath archived = work_dir /
            ("RaMAx.restart." + std::to_string(index) + ".log");
        if (std::filesystem::exists(archived)) continue;
        std::error_code error;
        std::filesystem::rename(current, archived, error);
        if (error) {
            throw std::runtime_error("Failed to archive previous RaMAx.log: " +
                                     error.message());
        }
        return;
    }
}

void clearNonReusableAlignmentState(const FilePath& work_dir) {
    constexpr std::array<const char*, 3> volatile_directories{
        RESULT_DIR, "mask_interval", "minipoa_tmp"};
    for (const char* name : volatile_directories) {
        const FilePath path = work_dir / name;
        std::error_code error;
        std::filesystem::remove_all(path, error);
        if (error) {
            throw std::runtime_error(
                "Failed to clear non-reusable restart state " + path.string() +
                ": " + error.message());
        }
    }
    spdlog::info(
        "Restart cache boundary applied: raw/clean/index retained; result/mask_interval/minipoa_tmp cleared");
    spdlog::info(
        "Alignment restarted from beginning; anchors, clusters, Blocks, DP state, and graph are not restored");
}

void configureOutputs(CommonArgs& args,
                      const std::vector<FilePath>& output_paths) {
    try {
        args.outputs = RaMAxOutput::validateOutputPaths(output_paths);
    } catch (const std::invalid_argument& error) {
        throw std::runtime_error(error.what());
    }
    args.output_paths = output_paths;
    args.output_path = output_paths.front();
    args.output_format = detectMultipleGenomeOutputFormat(args.output_path);
}

bool hasOutputFormat(const CommonArgs& args,
                     MultipleGenomeOutputFormat format) {
    return RaMAxOutput::hasFormat(args.outputs, format);
}

void validateHalAppendDependencyForOutputs(const CommonArgs& args) {
    const bool hal_output_requested =
        hasOutputFormat(args, MultipleGenomeOutputFormat::HAL);
    const auto executable =
        RaMAxDependencies::locateHalAppendCactusSubtreeExecutable();
    RaMAxDependencies::validateHalAppendCactusSubtree(
        executable, hal_output_requested);

    if (executable.empty()) {
        spdlog::warn(
            "[dependency-preflight] halAppendCactusSubtree was not found; "
            "continuing because no HAL output was requested");
    } else {
        spdlog::info(
            "[dependency-preflight] halAppendCactusSubtree={}",
            executable.string());
    }
}

std::string requiredMinipoaExecutable() {
    const auto executable =
        RaMesh::Alignment::locateMinipoaExecutable();
    if (executable.empty()) {
        throw std::runtime_error(
            "minipoa is required but was not found next to ramax or in PATH");
    }
    return executable.string();
}

}  // namespace


// ------------------------------------------------------------------
// 美化参数输出函数
// ------------------------------------------------------------------
inline void printRunConfiguration(const CommonArgs& args) {
    spdlog::info("");
    spdlog::info("============================================================");
    spdlog::info("                     RUN CONFIGURATION                     ");
    spdlog::info("============================================================");

    // Input/Output section
    spdlog::info("Input/Output:");
    spdlog::info("  Input file       : {}", args.input_path.string());
    spdlog::info("  Output files     : {}", args.outputs.size());
    for (const auto& output : args.outputs) {
        spdlog::info("    {:<3}           : {}",
            RaMAxOutput::formatName(output.format), output.path.string());
    }
    spdlog::info("  Work directory   : {}", args.work_dir_path.string());
    if (hasOutputFormat(args, MultipleGenomeOutputFormat::PAF)) {
        spdlog::info("  PAF mode         : {}", args.paf_mode);
    }
    if (hasOutputFormat(args, MultipleGenomeOutputFormat::GFA)) {
        spdlog::info("  GFA version      : {}", args.gfa_version);
        spdlog::info("  GFA profile      : {}", args.gfa_profile);
    }

    spdlog::info("");

    // Algorithm parameters section
    spdlog::info("Algorithm Parameters:");
    spdlog::info("  Chunk size            : {:L}", args.chunk_size);
    spdlog::info("  Overlap size          : {:L}", args.overlap_size);
    spdlog::info("  Min anchor length     : {}", args.min_anchor_length);
    spdlog::info("  Max anchor frequency  : {}", args.max_anchor_frequency);
    spdlog::info("  Search mode           : {}", SearchModeToString(args.search_mode));
    spdlog::info("  Accurate skip threshold: {}", args.accurate_skip_threshold);
    spdlog::info("  Allow MEM             : {}", args.allow_MEM ? "Enabled" : "Disabled");
    spdlog::info("  Fast build            : {}", args.fast_build ? "Enabled" : "Disabled");
    spdlog::info("  Reference coordinate sampling interval: {}",
                 args.sampling_interval);
    spdlog::info("  Suffix-array sampling rate            : {} (complete SA/ISA/LCP)",
                 args.sa_sampling_rate);
    spdlog::info("  Cluster Min span      : {}", args.min_span);
    spdlog::info("  Near distance         : {}", args.near_distance);
    spdlog::info("  Far distance          : {} (recorded only)", args.far_distance);
    spdlog::info("  Repeat masking        : {}", args.enable_repeat_masking ? "Enabled" : "Disabled");
    spdlog::info("  Tree root             ：{}", args.root_name);
    spdlog::info("  Ref genome name        : {}", args.ref_name.empty() ? "Not specified" : args.ref_name);
    spdlog::info("  Exact Block merge     : {}",
        args.merge_exact_contiguous_blocks ? "Enabled" : "Disabled");
    spdlog::info("  Query-gap merge max   : {}",
        args.merge_query_gap_max);
    spdlog::info("  Missing-species POA   : {}",
        args.realign_single_missing_species ? "Enabled" : "Disabled");
    if (args.realign_single_missing_species) {
        spdlog::info("  Missing-species span  : {}",
            args.species_mismatch_realign_max_span);
        spdlog::info("  Zero-gap merge span   : {}",
            args.species_mismatch_zero_gap_max_span);
    }
    spdlog::info("  Structural-break repair: {}",
        args.repair_structural_breaks ? "Enabled" : "Disabled");
    if (args.repair_structural_breaks) {
        spdlog::info("  Structural-break span: {}",
                     args.structural_break_max_span);
    }
    spdlog::info("  Short-Block repair     : {}",
        args.repair_short_blocks ? "Enabled" : "Disabled");
    spdlog::info("  Single round alignment: {}", args.one_round ? "Enabled" : "Disabled");
    spdlog::info("");


    // Performance section
    spdlog::info("Performance:");
    spdlog::info("  Thread count          : {}", args.thread_num);
    spdlog::info("  Restart mode          : {}", args.restart ? "Enabled" : "Disabled");


    spdlog::info("");

    // Output control section
    spdlog::info("Output Control:");
    spdlog::info("  Log level             : {}", args.log_level);


    spdlog::info("============================================================");
    spdlog::info("");
}

// ------------------------------------------------------------------
// CLI11 参数注册函数
// 用于配置和注册通用命令行参数（参考 RaMAx 主程序）
// ------------------------------------------------------------------
inline void setupCommonOptions(CLI::App* cmd, CommonArgs& args) {

    // ========================
    // CLI 帮助信息格式设置
    // ========================

    auto fmt = std::make_shared<CustomFormatter>(); // 自定义 CLI 帮助输出格式
    fmt->column_width(50);                          // 设置帮助信息列宽
    cmd->formatter(fmt);                            // 应用自定义格式化器

    // 版本信息参数（-v / --version）
    cmd->set_version_flag("-v,--version",
        std::string("RaMAx version ") + RAMAX_VERSION);

    // ========================
    // 输入 / 输出路径参数
    // ========================

    // 输入序列文件路径
    auto* input_opt = cmd->add_option("-i,--input", args.input_path,
        "Seqfile path; a Newick tree is required only for HAL output.")
        ->group("Input Files")                      // 帮助信息分组
        ->type_name("<path>")                       // 参数类型显示名称
        ->transform(trim_whitespace);               // 去除首尾空白字符

    // 输出结果路径
    cmd->add_option("-o,--output", args.output_paths,
        "Output alignment path; repeat -o for MAF, PAF, GFA, and HAL.")
        ->group("Output")
        ->type_name("<path>")
        ->transform(trim_whitespace);

    cmd->add_option("--paf-mode", args.paf_mode,
        "PAF pairing mode: connected or all (default: connected).")
        ->default_val("connected")
        ->capture_default_str()
        ->group("Output")
        ->type_name("<mode>")
        ->transform(CLI::CheckedTransformer(
            std::map<std::string, std::string>{
                {"connected", "connected"},
                {"all", "all"}
            },
            CLI::ignore_case));

    cmd->add_option("--gfa-version", args.gfa_version,
        "GFA path encoding: 1.0 uses P-lines; 1.1 uses W-lines (default: 1.1).")
        ->default_val("1.1")
        ->capture_default_str()
        ->group("Output")
        ->type_name("<version>")
        ->check(CLI::IsMember({"1.0", "1.1"}));

    cmd->add_option("--gfa-profile", args.gfa_profile,
        "GFA graph profile: exact preserves the audit graph; compact applies compact-v2-balanced (default: exact).")
        ->default_val("exact")
        ->capture_default_str()
        ->group("Output")
        ->type_name("<profile>")
        ->check(CLI::IsMember({"exact", "compact"}));

    // 工作目录路径（中间文件、索引缓存等）
    auto* workspace_opt = cmd->add_option("-w,--workdir", args.work_dir_path,
        "Working directory for intermediate files and logs.")
        ->group("Output")
        ->type_name("<path>")
        ->transform(trim_whitespace);

    // ========================
    // 基因组切片与锚点相关参数
    // ========================

    // 每个切片的长度（用于并行处理）
    cmd->add_option("--chunk_size", args.chunk_size,
        "Size of each chunk for parallel processing (default: 10000000).")
        ->default_val(10000000)                     // 默认值
        ->capture_default_str()                     // 在 --help 中显示默认值
        ->group("Software Parameters")
        ->check(CLI::Range(1000000,
            std::numeric_limits<int>::max()))       // 合法取值范围
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // HAL 文件中的根基因组名称
    cmd->add_option("--root", args.root_name,
        "HAL root genome name; valid only for HAL output.")
        ->group("Output")
        ->type_name("<string>")
        ->transform(trim_whitespace);

    // 参考基因组名称
    cmd->add_option("--ref", args.ref_name,
        "Ref genome name used in alignment")
        ->group("Software Parameters")
        ->type_name("<string>")
        ->transform(trim_whitespace);

    // 相邻切片之间的重叠长度
    cmd->add_option("--overlap_size", args.overlap_size,
        "Size of overlap between chunks (default: 100000).")
        ->default_val(0)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 锚点的最小长度
    cmd->add_option(
        "--min_anchor_length", args.min_anchor_length,
        "Minimum anchor length (default: 20).")
        ->default_val(20)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 锚点最大出现频率过滤阈值
    cmd->add_option(
        "--max_anchor_frequency", args.max_anchor_frequency,
        "Maximum anchor frequency filter (default: 50).")
        ->default_val(50)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // ========================
    // 搜索 / 算法参数
    // ========================

    // 锚点搜索模式（fast / middle / accurate）
    cmd->add_option("--search-mode", args.search_mode,
        "Anchor search mode: fast/middle/accurate (default: accurate).")
        ->default_val(ACCURATE_SEARCH)
        ->capture_default_str()
        ->group("Software Parameters")
        ->type_name("<mode>")
        ->transform(CLI::CheckedTransformer(
            std::map<std::string, SearchMode>{
                {"fast", FAST_SEARCH},
                {"middle", MIDDLE_SEARCH},
                {"accurate", ACCURATE_SEARCH}
            },
            CLI::ignore_case));                     // 忽略大小写

    cmd->add_option(
        "--accurate-skip-threshold", args.accurate_skip_threshold,
        "Skip accepted unique MUMs longer than this many bp in accurate mode; 0 disables (default: 10000).")
        ->default_val(5000)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 是否允许使用 MEM（Maximal Exact Match）
    cmd->add_flag("--allow-mem", args.allow_MEM,
        "Allow MEM (Maximal Exact Match) instead of only MUM.")
        ->group("Software Parameters");

    // 是否仅运行一轮比对流程
    cmd->add_flag("--one-round", args.one_round,
        "Only run one round for alignment.")
        ->group("Software Parameters");

    cmd->add_flag(
        "--optimize-blocks",
        "Explicitly enable the default Block optimizations (enabled by default).")
        ->group("Graph Optimization");

    cmd->add_flag(
        "--merge-blocks",
        args.merge_exact_contiguous_blocks,
        "Merge compatible neighboring Blocks.")
        ->group("Graph Optimization");

    cmd->add_option(
        "--merge-gap",
        args.merge_query_gap_max,
        "Maximum query gap allowed when merging Blocks (bp).")
        ->default_val(100)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->check(CLI::Range(0, 10000));

    auto* realign_missing_species_flag = cmd->add_flag(
        "--realign-missing",
        args.realign_single_missing_species,
        "Realign bounded windows with missing sequences.")
        ->group("Graph Optimization");

    cmd->add_option(
        "--realign-span",
        args.species_mismatch_realign_max_span,
        "Maximum span for missing-sequence realignment (bp).")
        ->default_val(3000)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->needs(realign_missing_species_flag)
        ->check(CLI::Range(1, 10000));

    cmd->add_option(
        "--zero-gap-span",
        args.species_mismatch_zero_gap_max_span,
        "Maximum span for zero-gap missing windows (bp).")
        ->default_val(200)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->check(CLI::Range(1, 3000));

    cmd->add_flag(
        "--repair-breaks",
        args.repair_structural_breaks,
        "Repair high-confidence structural discontinuities.")
        ->group("Graph Optimization");

    cmd->add_option(
        "--break-span",
        args.structural_break_max_span,
        "Maximum structural-break repair span (bp).")
        ->default_val(1000)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->check(CLI::Range(1, 1000));

    cmd->add_flag(
        "--merge-short-blocks",
        args.repair_short_blocks,
        "Try to merge short Blocks with banded KSW2.")
        ->group("Graph Optimization");

    // 使用慢但更精确的索引构建方式
    cmd->add_flag("--slow-build",
        "Use slow but more accurate index building method.")
        ->group("Software Parameters");

    // 参考序列索引采样间隔
    cmd->add_option("--sampling-interval",
        args.sampling_interval,
        "Reference coordinate-cache sampling interval (default: 32).")
        ->default_val(32)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);


    cmd->add_option("--sa-sampling-rate",
        args.sa_sampling_rate,
        "Suffix-array sampling rate; only the complete-array value 1 is "
        "supported by the hybrid suffix-array backend (default: 1).")
        ->default_val(1)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1, 1))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 构图或链式连接的最小跨度阈值
    cmd->add_option("--min-span", args.min_span,
        "Minimum span threshold for graph construction (default: 50).")
        ->default_val(65)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    cmd->add_option("--near-distance", args.near_distance,
        "Mash distance below which first-round queries use wfmash (default: 0.01).")
        ->default_val(0.01)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0.0, 1.0))
        ->type_name("<float>");

    cmd->add_option("--far-distance", args.far_distance,
        "Reserved distant-species Mash threshold; recorded only (default: 0.02).")
        ->default_val(0.02)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0.0, 1.0))
        ->type_name("<float>");

    // Repeat masking is currently disabled because it depends on external
    // windowmasker binaries that are no longer bundled with RaMAx.
    // cmd->add_flag("--mask-repeats", args.enable_repeat_masking,
    //     "Enable repeat sequence masking.")
    //     ->group("Software Parameters");

    // ========================
    // 性能相关参数
    // ========================

    // 并行线程数
    cmd->add_option("-t,--threads", args.thread_num,
        "Number of threads to use for parallel processing (default: system cores).")
        ->default_val(std::thread::hardware_concurrency())
        ->envname("RAMAx_THREADS")                  // 支持环境变量设置
        ->capture_default_str()
        ->group("Performance")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 复用预处理和锚点索引缓存；比对与构图始终重新开始
    auto* restart_flag = cmd->add_flag("--restart", args.restart,
        "Reuse raw/clean FASTA and anchor-index caches, then rerun alignment from the beginning; explicitly supplied options override saved values.")
        ->group("Performance");

    // ========================
    // 日志与输出控制参数
    // ========================

    // 日志级别设置
    cmd->add_option("--log-level", args.log_level,
        "Log level: debug/info/warn/error (default: info).")
        ->default_val("info")
        ->capture_default_str()
        ->group("Output Control")
        ->type_name("<level>")
        ->transform(CLI::CheckedTransformer(
            std::map<std::string, std::string>{
                {"debug", "debug"},
                {"info", "info"},
                {"warn", "warn"},
                {"error", "error"}
            },
            CLI::ignore_case));

    // 详细输出模式
    auto* verbose_flag = cmd->add_flag("--verbose", args.verbose,
        "Enable verbose output mode.")
        ->group("Output Control");

    // 静默模式（仅输出错误信息）
    auto* quiet_flag = cmd->add_flag("--quiet", args.quiet,
        "Enable quiet mode (only errors).")
        ->group("Output Control");

    // ========================
    // 参数依赖与互斥关系
    // ========================

    // --restart 必须依赖工作目录存在
    restart_flag->needs(workspace_opt);

    // 输入快照是 workdir 身份的一部分，restart 不允许替换 seqfile。
    restart_flag->excludes(input_opt);

    // verbose 与 quiet 互斥
    verbose_flag->excludes(quiet_flag);
    quiet_flag->excludes(verbose_flag);

}


// ------------------------------
// 工具函数：处理 --slow-build 标志
// ------------------------------
static inline void applySlowBuildFlag(const CLI::App& app, CommonArgs& common_args) {
    // 若用户指定 --slow-build，则关闭 fast_build
    if (app.count("--slow-build")) {
        common_args.fast_build = false;
    }
}

static inline void applyGraphOptimizationOptions(
    const CLI::App& app, CommonArgs& args) {
    if (app.count("--optimize-blocks") != 0) {
        args.merge_exact_contiguous_blocks = true;
        args.realign_single_missing_species = true;
        args.repair_structural_breaks = true;
        args.repair_short_blocks = true;
    }

    const bool merge_enabled = args.merge_exact_contiguous_blocks;
    const bool realign_enabled = args.realign_single_missing_species;
    const bool break_enabled = args.repair_structural_breaks;
    if (realign_enabled && !merge_enabled) {
        throw CLI::ValidationError(
            "--realign-missing requires --merge-blocks");
    }
    if (app.count("--merge-gap") != 0 && !merge_enabled) {
        throw CLI::ValidationError("--merge-gap requires --merge-blocks");
    }
    if ((app.count("--realign-span") != 0 ||
         app.count("--zero-gap-span") != 0) &&
        !realign_enabled) {
        throw CLI::ValidationError(
            "--realign-span and --zero-gap-span require --realign-missing");
    }
    if (app.count("--break-span") != 0 && !break_enabled) {
        throw CLI::ValidationError("--break-span requires --repair-breaks");
    }
}

// ------------------------------
// 工具函数：根据 quiet/verbose/log_level 设置 spdlog 日志级别
// ------------------------------
static inline void configureLogLevel(const CommonArgs& common_args) {
    // quiet 优先：只输出 error
    if (common_args.quiet) {
        spdlog::set_level(spdlog::level::err);
        return;
    }

    // verbose 次之：输出 debug
    if (common_args.verbose) {
        spdlog::set_level(spdlog::level::debug);
        return;
    }

    // 否则按照 log_level 字符串设置
    if (common_args.log_level == "debug") {
        spdlog::set_level(spdlog::level::debug);
    } else if (common_args.log_level == "info") {
        spdlog::set_level(spdlog::level::info);
    } else if (common_args.log_level == "warn") {
        spdlog::set_level(spdlog::level::warn);
    } else if (common_args.log_level == "error") {
        spdlog::set_level(spdlog::level::err);
    }
}

static inline void logBuildMode() {
#ifdef _DEBUG_
    spdlog::info("Build mode: Debug");
#else
    spdlog::info("Build mode: Release");
#endif
}

// ------------------------------
// 工具函数：确保工作目录存在且合法；必要时创建
//（用于 restart / normal 两种模式的目录检查）
// ------------------------------
static inline void ensureWorkDirValidOrCreate(const std::filesystem::path& work_dir_path) {
    // 工作目录存在：必须是目录
    if (std::filesystem::exists(work_dir_path)) {
        if (!std::filesystem::is_directory(work_dir_path)) {
            throw CLI::ValidationError("Work directory is not valid: " + work_dir_path.string());
        }
    }
    // 工作目录不存在：创建
    else {
        std::filesystem::create_directories(work_dir_path);
    }
}

// ------------------------------
// 模式 1：重启模式（--restart）
// - workdir 必须提供
// - 读取 workdir/CONFIG_FILE 反序列化 common_args
// ------------------------------
static int runRestartMode(CommonArgs& common_args,
                          const RestartOverrides& overrides,
                          const CLI::App& app) {
    // restart 模式下 workdir 必须提供
    if (common_args.work_dir_path.empty()) {
        throw CLI::RequiredError("In restart mode, --workdir (-w) is required.");
    }

    const FilePath requested_work_dir = common_args.work_dir_path;
    if (!std::filesystem::is_directory(requested_work_dir)) {
        throw CLI::ValidationError(
            "Restart work directory does not exist: " +
            requested_work_dir.string());
    }

    archivePreviousLog(requested_work_dir);
    setupLoggerWithFile(requested_work_dir);
    configureLogLevel(common_args);
    spdlog::info("RaMAx version {}", RAMAX_VERSION);
    logBuildMode();
    spdlog::info("Restart mode enabled.");

    // 加载之前保存的参数配置文件
    const FilePath config_path = requested_work_dir / CONFIG_FILE;
    const uint32_t schema = detectConfigSchema(config_path);
    CommonArgs loaded;
    if (schema == CONFIG_SCHEMA_VERSION) {
        try {
            std::ifstream input(config_path);
            cereal::JSONInputArchive archive(input);
            archive(cereal::make_nvp("common_args", loaded));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Restart configuration is incompatible with RaMAx " +
                std::string(RAMAX_VERSION) + ": " + error.what());
        }
    } else if (schema == 5) {
        LegacyCommonArgsV5 legacy;
        try {
            std::ifstream input(config_path);
            cereal::JSONInputArchive archive(input);
            archive(cereal::make_nvp("common_args", legacy));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Invalid schema-5 restart configuration: " +
                std::string(error.what()));
        }
        loaded = convertSchema5Config(legacy);
        spdlog::warn(
            "Loaded schema-5 workdir: --sa-sampling-rate=1 was added "
            "during migration");
    } else if (schema == 4) {
        LegacyCommonArgsV4 legacy;
        try {
            std::ifstream input(config_path);
            cereal::JSONInputArchive archive(input);
            archive(cereal::make_nvp("common_args", legacy));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Invalid schema-4 restart configuration: " +
                std::string(error.what()));
        }
        loaded = convertSchema4Config(legacy);
        spdlog::warn(
            "Loaded schema-4 workdir: --gfa-profile=exact was added during migration");
    } else if (schema == 3) {
        LegacyCommonArgsV3 legacy;
        try {
            std::ifstream input(config_path);
            cereal::JSONInputArchive archive(input);
            archive(cereal::make_nvp("common_args", legacy));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Invalid schema-3 restart configuration: " +
                std::string(error.what()));
        }
        loaded = convertSchema3Config(legacy);
        spdlog::warn(
            "Loaded schema-3 workdir: --gfa-version=1.1 was added during migration");
    } else if (schema == 2) {
        LegacyCommonArgsV2 legacy;
        try {
            std::ifstream input(config_path);
            cereal::JSONInputArchive archive(input);
            archive(cereal::make_nvp("common_args", legacy));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Invalid schema-2 restart configuration: " +
                std::string(error.what()));
        }
        loaded = convertSchema2Config(legacy);
        spdlog::warn(
            "Loaded schema-2 workdir: --near-distance=0.01 and --far-distance=0.02 were added during migration");
    } else if (schema == 1) {
        LegacyCommonArgsV1 legacy;
        try {
            std::ifstream input(config_path);
            cereal::JSONInputArchive archive(input);
            archive(cereal::make_nvp("common_args", legacy));
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Invalid schema-1 restart configuration: " +
                std::string(error.what()));
        }
        loaded = convertLegacyConfig(legacy);
        std::vector<FilePath> legacy_outputs{loaded.output_path};
        const FilePath outputs_path = requested_work_dir / OUTPUTS_FILE;
        if (std::filesystem::exists(outputs_path)) {
            OutputManifest manifest;
            try {
                std::ifstream input(outputs_path);
                cereal::JSONInputArchive archive(input);
                archive(cereal::make_nvp("outputs", manifest));
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "Invalid legacy output manifest: " +
                    std::string(error.what()));
            }
            if (manifest.schema_version != OUTPUTS_SCHEMA_VERSION ||
                manifest.output_paths.empty()) {
                throw std::runtime_error("Invalid legacy output manifest schema");
            }
            legacy_outputs = std::move(manifest.output_paths);
        }
        loaded.output_paths = std::move(legacy_outputs);

        const FilePath threshold_path =
            requested_work_dir / "accurate_skip_threshold.txt";
        if (std::filesystem::exists(threshold_path)) {
            std::ifstream input(threshold_path);
            uint64_t value = 0;
            if (!(input >> value) || value > std::numeric_limits<uint_t>::max()) {
                throw std::runtime_error("Invalid legacy accurate skip threshold");
            }
            input >> std::ws;
            if (!input.eof()) {
                throw std::runtime_error("Invalid legacy accurate skip threshold");
            }
            loaded.accurate_skip_threshold = static_cast<uint_t>(value);
        }
        const FilePath paf_mode_path = requested_work_dir / "paf_mode.txt";
        if (std::filesystem::exists(paf_mode_path)) {
            std::ifstream input(paf_mode_path);
            if (!(input >> loaded.paf_mode) ||
                (loaded.paf_mode != "connected" && loaded.paf_mode != "all")) {
                throw std::runtime_error("Invalid legacy PAF mode");
            }
        }
        spdlog::warn(
            "Loaded schema-1 workdir: missing fields use RaMAx 1.0.6 restart defaults; existing caches will be trusted once and migrated");
    } else {
        throw std::runtime_error(
            "Unsupported restart schema version " + std::to_string(schema));
    }

    loaded.work_dir_path = requested_work_dir;
    loaded.restart = true;
    loaded.pending_config_update = true;
    applyRestartOverrides(loaded, overrides);
    applyGraphOptimizationOptions(app, loaded);
    configureOutputs(loaded, loaded.output_paths);
    validateHalAppendDependencyForOutputs(loaded);

    if (loaded.overlap_size >= loaded.chunk_size) {
        throw std::runtime_error("Overlap size must be less than chunk size.");
    }
    validateDistanceThresholds(loaded);
    validateSaSamplingRate(loaded);
    validateGfaVersion(loaded);
    validateGfaProfile(loaded);
    if (loaded.paf_mode_explicit &&
        !hasOutputFormat(loaded, MultipleGenomeOutputFormat::PAF)) {
        throw std::runtime_error("--paf-mode is only valid with .paf output");
    }
    if (loaded.gfa_version_explicit &&
        !hasOutputFormat(loaded, MultipleGenomeOutputFormat::GFA)) {
        throw std::runtime_error(
            "--gfa-version is only valid with .gfa output");
    }
    if (loaded.gfa_profile_explicit &&
        !hasOutputFormat(loaded, MultipleGenomeOutputFormat::GFA)) {
        throw std::runtime_error(
            "--gfa-profile is only valid with .gfa output");
    }
    if (overrides.has("--root") &&
        !hasOutputFormat(loaded, MultipleGenomeOutputFormat::HAL)) {
        throw std::runtime_error("--root is only supported for HAL output");
    }
    if (!hasOutputFormat(loaded, MultipleGenomeOutputFormat::HAL)) {
        loaded.root_name.clear();
    }

    common_args = std::move(loaded);
    configureLogLevel(common_args);
    spdlog::info("Effective restart configuration loaded from {}",
                 config_path.string());
    spdlog::info("Explicit restart overrides: {}", overrides.specified.size());
    return 0;
}

// ------------------------------
// 模式 2：正常运行模式
// - 检查必需参数 input/output/workdir
// - 非 debug 下：工作目录必须为空（否则报错）
// - 保存参数到 workdir/CONFIG_FILE（用于 --restart）
// ------------------------------
static int runNormalMode(CommonArgs& common_args, const CLI::App& app) {
    // 检查必要参数
    if (common_args.input_path.empty())
        throw CLI::RequiredError("Missing required option: --input (-i)");
    if (common_args.output_paths.empty())
        throw CLI::RequiredError("Missing required option: --output (-o)");
    if (common_args.work_dir_path.empty())
        throw CLI::RequiredError("Missing required option: --workdir (-w)");

    applyGraphOptimizationOptions(app, common_args);
    configureOutputs(common_args, common_args.output_paths);
    validateHalAppendDependencyForOutputs(common_args);

    if (!hasOutputFormat(common_args, MultipleGenomeOutputFormat::PAF) &&
        common_args.paf_mode_explicit) {
        throw std::runtime_error(
            "--paf-mode is only valid with .paf output");
    }
    if (!hasOutputFormat(common_args, MultipleGenomeOutputFormat::GFA) &&
        common_args.gfa_version_explicit) {
        throw std::runtime_error(
            "--gfa-version is only valid with .gfa output");
    }
    if (!hasOutputFormat(common_args, MultipleGenomeOutputFormat::GFA) &&
        common_args.gfa_profile_explicit) {
        throw std::runtime_error(
            "--gfa-profile is only valid with .gfa output");
    }
    validateDistanceThresholds(common_args);
    validateSaSamplingRate(common_args);
    validateGfaVersion(common_args);
    validateGfaProfile(common_args);

#ifndef _DEBUG_
    // 非调试模式：确保工作目录为空且合法
    if (std::filesystem::exists(common_args.work_dir_path)) {
        if (!std::filesystem::is_directory(common_args.work_dir_path)) {
            throw CLI::ValidationError("Work directory is not valid: " + common_args.work_dir_path.string());
        }
        if (!std::filesystem::is_empty(common_args.work_dir_path)) {
            throw CLI::ValidationError("Work directory is not empty: " + common_args.work_dir_path.string());
        }
    } else {
        std::filesystem::create_directories(common_args.work_dir_path);
    }
#else
    // Debug 模式：若不存在则创建（保持原逻辑结构，不新增约束）
    ensureWorkDirValidOrCreate(common_args.work_dir_path);
#endif

    // 初始化日志器（输出到文件）
    setupLoggerWithFile(common_args.work_dir_path);
    configureLogLevel(common_args);
    spdlog::info("RaMAx version {}", RAMAX_VERSION);
    logBuildMode();
    spdlog::info("Multiple genome alignment mode enabled.");

    // 检查 chunk 与 overlap 的合法性
    if (common_args.overlap_size >= common_args.chunk_size) {
        throw std::runtime_error("Overlap size must be less than chunk size.");
    }
    common_args.schema_version = CONFIG_SCHEMA_VERSION;
    saveEffectiveConfig(common_args);
    const FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
    spdlog::info("Configuration saved to {}", config_path.string());

    return 0;
}

// ------------------------------
// 运行前的配置阶段：
// - 根据 restart 标志选择模式
// - 捕获 runtime_error 并统一打印提示
// ------------------------------
static int prepareRun(CommonArgs& common_args, const CLI::App& app,
                      const RestartOverrides& overrides) {
    try {
        // 模式 1：重启模式（--restart）
        if (common_args.restart) {
            return runRestartMode(common_args, overrides, app);
        }

        // 模式 2：正常运行模式
        return runNormalMode(common_args, app);
    }
    catch (const std::runtime_error& e) {
        spdlog::error("{}", e.what());
        spdlog::error("Use --help for usage information.");
        spdlog::error("Exiting with error code 1.");
        return 1;
    }
}

// ------------------------------
// 输入验证阶段：解析 seqfile 与校验输入基因组路径（本地/URL）
// ------------------------------
static int inputValidationPhase(
    CommonArgs& common_args,
    NewickParser& newick_tree,
    SpeciesPathMap& species_path_map
) {
    if (!hasOutputFormat(common_args, MultipleGenomeOutputFormat::HAL) &&
        !common_args.root_name.empty()) {
        throw std::runtime_error("--root is only supported for HAL output");
    }

    // The input identity is independent of a HAL subtree selection, so parse
    // the complete seqfile before applying --root.
    SpeciesPathMap all_species;
    NewickParser complete_tree;
    const bool complete_has_tree = parseSeqfile(
        common_args.input_path, complete_tree, all_species, "");
    validateOrCreateInputManifest(common_args, all_species);

    // HAL requires a species tree. Other output formats may use mappings only.
    std::string root = common_args.root_name;
    const bool has_tree =
        parseSeqfile(common_args.input_path, newick_tree, species_path_map, root);
    if (has_tree != complete_has_tree) {
        throw std::runtime_error("Seqfile tree parsing changed during validation");
    }
    if (hasOutputFormat(common_args, MultipleGenomeOutputFormat::HAL) &&
        !has_tree) {
        throw std::runtime_error(
            "HAL output requires a Newick tree as the first seqfile record");
    }

    // 逐个校验输入基因组路径是否合法（URL 可达 / 本地文件存在）
    for (const auto& [species, path] : species_path_map) {
        if (isUrl(path.string())) {
            if (common_args.restart &&
                hasReusableRawSnapshot(common_args, species, path)) {
                spdlog::info(
                    "Using cached raw snapshot without rechecking URL: {}",
                    species);
            } else {
                verifyUrlReachable(path.string());
            }
        } else {
            verifyLocalFile(path);
        }

        spdlog::info("Input genome: {} (size: {})",
                     species,
                     getReadableFileSize(path));
    }

    if (common_args.pending_config_update) {
        saveEffectiveConfig(common_args);
        common_args.pending_config_update = false;
        spdlog::info("Restart configuration migrated/updated to schema {}",
                     CONFIG_SCHEMA_VERSION);
    }

    return 0;
}

// ------------------------------
// 数据预处理阶段：
// - copyRawData（拷贝/下载）
// - cleanRawDataset（清洗）
// - 可选 repeat masking 生成 interval
// - 创建 SeqPro managers（带/不带 interval）
// ------------------------------
static int preprocessingPhase(
    CommonArgs& common_args,
    SpeciesPathMap& species_path_map,
    std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    SoftMask::PathMap& softmask_path_map,
    SeqPro::Length& reference_min_seq_length
) {
    // 拷贝或下载原始文件（并行执行）
    RaMAxCache::StageStats raw_stats;
    copyRawData(common_args.work_dir_path, species_path_map,
                common_args.thread_num, &raw_stats,
                common_args.trust_legacy_cache);

    // interval 文件映射（重复遮蔽用）
    std::map<SpeciesName, FilePath> interval_files_map;

    // 参考序列最短长度（用于后续 sampling_interval 截断）
    reference_min_seq_length = std::numeric_limits<SeqPro::Length>::max();

    const bool align_softmasked_regions =
        std::getenv("RAMAX_ALIGN_SOFTMASK") != nullptr;
    // HAL always records the original lowercase runs for lossless export.
    // The experimental alignment switch additionally loads those runs as
    // seed-mask intervals while retaining the complete original coordinate
    // space for extension and export. Cache accounting remains identical for
    // either preprocessing route.
    RaMAxCache::StageStats clean_stats;
    if (hasOutputFormat(common_args, MultipleGenomeOutputFormat::HAL) ||
        align_softmasked_regions) {
        cleanRawDatasetWithSoftMaskIndex(
            common_args.work_dir_path, species_path_map,
            softmask_path_map, common_args.thread_num,
            &clean_stats, common_args.trust_legacy_cache);
    } else {
        cleanRawDataset(
            common_args.work_dir_path, species_path_map,
            common_args.thread_num, &clean_stats,
            common_args.trust_legacy_cache);
    }
    spdlog::info(
        "[cache-summary] preprocessing raw(reused/rebuilt)={}/{} clean(reused/rebuilt)={}/{}",
        raw_stats.reused, raw_stats.rebuilt,
        clean_stats.reused, clean_stats.rebuilt);

    if (common_args.enable_repeat_masking) {
        spdlog::warn("Repeat masking is currently disabled; continuing without repeat masking.");
        common_args.enable_repeat_masking = false;
    }

    // 若启用重复遮蔽：生成 interval 并创建 MaskedSequenceManager(带 interval)
    if (common_args.enable_repeat_masking) {
        spdlog::info("Repeat masking enabled. Generating interval files based on raw files...");

        // 生成 interval 文件
        interval_files_map = repeatSeqMasking(
            common_args.work_dir_path, species_path_map, common_args.thread_num);

        if (interval_files_map.empty()) {
            spdlog::error("No interval files were generated. Please check the input FASTA files and ensure they are valid.");
            return 1;
        }

        spdlog::info("Interval files generated successfully.");
        spdlog::info("Creating SeqPro managers with repeat masking support...");

        for (const auto& [species_name, cleaned_fasta_path] : species_path_map) {
            if (interval_files_map.contains(species_name)) {
                try {
                    auto original_manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
                    auto manager = std::make_unique<SeqPro::MaskedSequenceManager>(
                        std::move(original_manager),
                        interval_files_map[species_name]
                    );

                    spdlog::info("[{}] SeqPro Manager created with repeat masking: {}",
                                 species_name, cleaned_fasta_path.string());

                    // 记录序列统计信息
                    SequenceUtils::recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);

                    auto shared_manager = std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
                    seqpro_managers[species_name] = std::move(shared_manager);
                }
                catch (const std::exception& e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name, e.what());
                    return 1;
                }
            }
            else {
                try {
                    auto manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);

                    spdlog::info("[{}] SeqPro Manager created without repeat masking: {}",
                                 species_name, cleaned_fasta_path.string());

                    // 记录序列统计信息
                    SequenceUtils::recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);

                    auto shared_manager = std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
                    seqpro_managers[species_name] = std::move(shared_manager);
                }
                catch (const std::exception& e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name, e.what());
                    return 1;
                }
            }
        }
    }
    // 不启用重复遮蔽：基于清洗后的文件创建 MaskedSequenceManager（无 interval）
    else {
        spdlog::info("Repeat masking disabled. Creating standard SeqPro managers...");

        for (const auto& [species_name, cleaned_fasta_path] : species_path_map) {
            try {
                auto original_manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
                auto manager = std::make_unique<SeqPro::MaskedSequenceManager>(std::move(original_manager));
                if (align_softmasked_regions) {
                    const auto softmask_it =
                        softmask_path_map.find(species_name);
                    if (softmask_it == softmask_path_map.end()) {
                        throw std::runtime_error(
                            "Missing soft-mask index for " + species_name);
                    }
                    const SoftMask::Index softmask_index(
                        softmask_it->second);
                    uint64_t masked_bases = 0;
                    for (const auto& sequence_name :
                         manager->getSequenceNames()) {
                        const auto indexed_intervals =
                            softmask_index.intervals(sequence_name);
                        std::vector<SeqPro::MaskInterval> intervals;
                        intervals.reserve(indexed_intervals.size());
                        for (const auto& [start, end] :
                             indexed_intervals) {
                            intervals.emplace_back(start, end);
                            masked_bases += end - start;
                        }
                        manager->addMaskIntervals(
                            sequence_name, intervals);
                    }
                    manager->finalizeMaskIntervals();
                    spdlog::info(
                        "[{}] Loaded {} soft-masked bases for "
                        "alignment seeding",
                        species_name, masked_bases);
                }

                spdlog::info("[{}] SeqPro Manager created: {}", species_name, cleaned_fasta_path.string());

                // 记录序列统计信息
                SequenceUtils::recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);

                auto shared_manager = std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
                seqpro_managers[species_name] = std::move(shared_manager);
            }
            catch (const std::exception& e) {
                spdlog::error("[{}] Error creating SeqPro manager: {}", species_name, e.what());
                return 1;
            }
        }
    }

    return 0;
}

// ------------------------------
// 清理 SeqPro managers 的遮蔽区间（导出前）
// ------------------------------
static void clearAllMaskedRegions(
    const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers
) {
    for (const auto& [species_name, seq_mgr_variant] : seqpro_managers) {
        std::visit([&](const auto& seq_mgr) {
            seq_mgr->clearMaskedRegions();
        }, *seq_mgr_variant);
    }
}

// ------------------------------
// 运行 star alignment，并返回构建好的图对象
// ------------------------------
static std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> runStarAlignment(
    const CommonArgs& common_args,
    SpeciesPathMap& species_path_map,
    std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    SeqPro::Length reference_min_seq_length,
    double& align_seconds_out
) {
    const std::string minipoa_executable = requiredMinipoaExecutable();
    RaMesh::Alignment::ExternalMsaRunner::instance()
        .configureScratchDirectory(
            common_args.work_dir_path / "minipoa_tmp");

    // 初始化比对器
    MultipleRareAligner mra(
        common_args.work_dir_path,
        species_path_map,
        common_args.thread_num,
        common_args.chunk_size,
        common_args.overlap_size,
        common_args.min_anchor_length,
        common_args.max_anchor_frequency,
        common_args.accurate_skip_threshold,
        common_args.trust_legacy_cache
    );
    mra.allow_mem = common_args.allow_MEM;
    mra.sa_sampling_rate = common_args.sa_sampling_rate;

    mra.merge_exact_contiguous_blocks_enabled =
        common_args.merge_exact_contiguous_blocks;
    mra.merge_query_gap_max = common_args.merge_query_gap_max;

    // 计时：star alignment 总耗时
    mra.realign_single_missing_species_enabled =
        common_args.realign_single_missing_species;
    mra.species_mismatch_realign_max_span =
        common_args.species_mismatch_realign_max_span;
    mra.species_mismatch_zero_gap_max_span =
        common_args.species_mismatch_zero_gap_max_span;
    mra.species_mismatch_msa_executable = minipoa_executable;
    mra.structural_break_repair_options.enabled =
        common_args.repair_structural_breaks;
    mra.structural_break_repair_options.maximum_span =
        common_args.structural_break_max_span;
    mra.structural_break_repair_options.parallel_threads =
        common_args.thread_num;
    mra.structural_break_repair_options.msa_executable = minipoa_executable;
    mra.short_block_repair_options.enabled =
        common_args.repair_short_blocks;
    mra.short_block_repair_options.maximum_query_gap =
        common_args.merge_query_gap_max;
    mra.short_block_repair_options.parallel_threads =
        common_args.thread_num;
    mra.near_distance_threshold = common_args.near_distance;
    mra.far_distance_threshold = common_args.far_distance;
    auto t_start_align = std::chrono::steady_clock::now();

    // 初始化采样间隔：确保不超过 reference_min_seq_length（避免越界/无效采样）
    auto sampling_interval = std::min(
        static_cast<SeqPro::Length>(common_args.sampling_interval),
        reference_min_seq_length
    );

    // 执行 star alignment，构建多基因组图
    std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> graph =
        mra.starAlignment(
            seqpro_managers,
            common_args.ref_name,
            common_args.one_round,
            common_args.fast_build,
            sampling_interval,
            common_args.min_span
        );

    auto t_end_align = std::chrono::steady_clock::now();
    std::chrono::duration<double> align_time = t_end_align - t_start_align;
    align_seconds_out = align_time.count();

    return graph;
}

// ------------------------------
// 导出结果（MAF / PAF / GFA / HAL）
// ------------------------------
struct ExportAttempt {
    RaMAxOutput::OutputSpec output;
    bool success{false};
    double elapsed_seconds{0.0};
    std::string error;
};

static bool exportResults(
    const CommonArgs& common_args,
    const NewickParser& newick_tree,
    std::map<SpeciesName, SeqPro::SharedManagerVariant>& seqpro_managers,
    const SoftMask::PathMap& softmask_path_map,
    RaMesh::RaMeshMultiGenomeGraph* graph
) {
    // 清理遮蔽区间（导出前）
    clearAllMaskedRegions(seqpro_managers);

  // Both direct MAF export and HAL construction reconstruct multiway
  // columns through mergeAlignmentByRef().  Configure the shared repair
  // once so the two output paths make the same local insertion decisions.
  configureCrossAnchorInsertionRepair(
      requiredMinipoaExecutable(),
      common_args.species_mismatch_realign_max_span);
  graph->materializeSecondaryAlignments();
  if (common_args.merge_exact_contiguous_blocks &&
      !common_args.ref_name.empty()) {
    const size_t eliminated_boundaries = graph->mergeExactContiguousBlocks(
        common_args.ref_name, 1000000, common_args.merge_query_gap_max);
    spdlog::info("[secondary-alignments] post-materialization "
                 "eliminated_boundaries={} max_query_gap={}",
                 eliminated_boundaries, common_args.merge_query_gap_max);
  }

    std::vector<ExportAttempt> attempts;
    attempts.reserve(common_args.outputs.size());
    bool all_succeeded = true;

    for (const auto& output : common_args.outputs) {
        const auto started = std::chrono::steady_clock::now();
        ExportAttempt attempt;
        attempt.output = output;
        try {
            switch (output.format) {
            case MultipleGenomeOutputFormat::MAF:
                spdlog::info("Exporting MAF to {}...", output.path.string());
                graph->exportToMaf(
                    output.path, seqpro_managers, false);
                break;

            case MultipleGenomeOutputFormat::PAF: {
                spdlog::info(
                    "Exporting PAF to {} (mode={})...",
                    output.path.string(), common_args.paf_mode);
                RaMesh::Paf::PafExportOptions options;
                options.mode = common_args.paf_mode == "all"
                    ? RaMesh::Paf::Mode::ALL
                    : RaMesh::Paf::Mode::CONNECTED;
                options.only_primary = true;
                graph->exportToPaf(output.path, seqpro_managers, options);
                break;
            }

            case MultipleGenomeOutputFormat::GFA: {
                spdlog::info(
                    "Exporting GFA {} profile={} to {}...",
                    common_args.gfa_version, common_args.gfa_profile,
                    output.path.string());
                RaMesh::Gfa::GfaExportOptions options;
                options.version = RaMesh::Gfa::parseVersion(
                    common_args.gfa_version);
                options.profile = RaMesh::Gfa::parseProfile(
                    common_args.gfa_profile);
                options.only_primary = true;
                options.threads = common_args.thread_num;
                options.work_dir = common_args.work_dir_path;
                graph->exportToGfa(output.path, seqpro_managers, options);
                break;
            }

            case MultipleGenomeOutputFormat::HAL:
                spdlog::info("Exporting HAL to {}...", output.path.string());
                graph->exportToHal(
                    output.path,
                    seqpro_managers,
                    newick_tree,
                    common_args.root_name,
                    static_cast<int>(common_args.thread_num),
                    softmask_path_map);
                break;

            case MultipleGenomeOutputFormat::UNKNOWN:
                throw std::runtime_error(
                    "Unsupported output format for multiple genome alignment");
            }
            attempt.success = true;
        } catch (const std::exception& error) {
            attempt.error = error.what();
            all_succeeded = false;
            spdlog::error(
                "{} export failed for {}: {}",
                RaMAxOutput::formatName(output.format),
                output.path.string(), attempt.error);
        }
        attempt.elapsed_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        attempts.push_back(std::move(attempt));
    }
    logCrossAnchorInsertionRepairStats();
    RaMesh::Alignment::ExternalMsaRunner::instance().logSummary();

    spdlog::info("Output export summary:");
    for (const auto& attempt : attempts) {
        if (attempt.success) {
            spdlog::info(
                "  {} SUCCESS path={} elapsed_seconds={:.3f}",
                RaMAxOutput::formatName(attempt.output.format),
                attempt.output.path.string(), attempt.elapsed_seconds);
        } else {
            spdlog::error(
                "  {} FAILED path={} elapsed_seconds={:.3f} error={}",
                RaMAxOutput::formatName(attempt.output.format),
                attempt.output.path.string(), attempt.elapsed_seconds,
                attempt.error);
        }
    }
    return all_succeeded;
}

// ------------------------------
// 主流程：输入校验 -> 预处理 -> 比对 -> 导出
// 保持原有异常处理策略：捕获 std::runtime_error 并返回 1
// ------------------------------
static int runMainPipeline(CommonArgs& common_args, int argc, char** argv) {
    // 打印运行配置
    printRunConfiguration(common_args);

    // 打印执行命令行
    spdlog::info("Executed command: {}", getCommandLine(argc, argv));
    spdlog::info("");
    spdlog::info("============================================================");
    spdlog::info("                      INPUT VALIDATION                     ");
    spdlog::info("============================================================");

    SpeciesPathMap species_path_map;
    NewickParser newick_tree;

    try {
        // ------------------------------
        // 输入校验阶段
        // ------------------------------
        if (inputValidationPhase(common_args, newick_tree, species_path_map) != 0) {
            return 1;
        }

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                    DATA PREPROCESSING                     ");
        spdlog::info("============================================================");

        // ------------------------------
        // 数据预处理阶段
        // ------------------------------
        std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers;
        SoftMask::PathMap softmask_path_map;
        SeqPro::Length reference_min_seq_length = std::numeric_limits<SeqPro::Length>::max();

        if (preprocessingPhase(common_args, species_path_map, seqpro_managers,
                               softmask_path_map, reference_min_seq_length) != 0) {
            return 1;
        }

        if (common_args.restart) {
            clearNonReusableAlignmentState(common_args.work_dir_path);
        }

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                    STAR ALIGNMENT                         ");
        spdlog::info("============================================================");

        // ------------------------------
        // Star alignment 阶段
        // ------------------------------
        double align_seconds = 0.0;
        std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> graph =
            runStarAlignment(common_args, species_path_map,
                             seqpro_managers, reference_min_seq_length, align_seconds);

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                       COMPLETION                          ");
        spdlog::info("============================================================");
        spdlog::info("Star alignment completed in {:.3f} seconds.", align_seconds);

        // ------------------------------
        // 导出阶段
        // ------------------------------
        const bool exports_succeeded = exportResults(
            common_args, newick_tree, seqpro_managers,
            softmask_path_map, graph.get());
        if (!exports_succeeded) {
            spdlog::error(
                "One or more output formats failed; preserving work directory: {}",
                common_args.work_dir_path.string());
            return 1;
        }

        // Keep the work directory because it now contains the requested Mash
        // routing table, samtools-generated FAI provenance, and final
        // per-query mapping/alignment PAF files.
        spdlog::info("Work directory preserved with routing artifacts: {}",
                     common_args.work_dir_path.string());

        // ------------------------------
        // 退出
        // ------------------------------
        spdlog::info("RaMAx execution completed successfully!");
        return 0;
    }
    catch (const std::runtime_error& e) {
        spdlog::error("{}", e.what());
        spdlog::error("Use --help for usage information.");
        spdlog::error("Exiting with error code 1.");
        return 1;
    }
}

// ------------------------------------------------------------------
// main：仅保留“初始化/解析/模式分发/主流程调用”的骨架
// ------------------------------------------------------------------
int main(int argc, char** argv) {
    // 初始化异步日志线程池（spdlog）
    // 日志缓冲区容量 8192，单线程写入（异步落盘）
    spdlog::init_thread_pool(8192, 1);

    // 初始化 CLI 应用
    CLI::App app{ "RaMAx: A High-performance Genome Alignment Tool" };

    // 存储用户输入参数
    CommonArgs common_args;

    // 注册通用参数
    setupCommonOptions(&app, common_args);

    // 开始解析命令行参数
    CLI11_PARSE(app, argc, argv);
    common_args.paf_mode_explicit = app.count("--paf-mode") != 0;
    common_args.gfa_version_explicit = app.count("--gfa-version") != 0;
    common_args.gfa_profile_explicit = app.count("--gfa-profile") != 0;
    applySlowBuildFlag(app, common_args);
    const RestartOverrides restart_overrides =
        captureRestartOverrides(app, common_args);

    // Configure early console diagnostics; restart applies the saved/effective
    // log policy after loading its configuration.
    configureLogLevel(common_args);

    // Resolve the three unconditional external dependencies before creating or
    // mutating the work directory. CLI11 has already handled --help and
    // --version, so those informational commands remain dependency-free. The
    // HAL-only helper is checked after the effective output list is known.
    try {
        const auto dependencies =
            RaMAxDependencies::requireUnconditionalStartupDependencies();
        spdlog::info(
            "[dependency-preflight] minipoa={}",
            dependencies.minipoa.string());
        spdlog::info(
            "[dependency-preflight] wfmash={}",
            dependencies.wfmash.string());
        spdlog::info(
            "[dependency-preflight] mash={}",
            dependencies.mash.string());
        spdlog::info(
            "[dependency-preflight] all unconditional required executables "
            "are available");
    } catch (const std::exception& error) {
        spdlog::critical("{}", error.what());
        return 1;
    }

    // 运行前准备：根据 restart 与否进行目录/参数/配置文件处理
    if (prepareRun(common_args, app, restart_overrides) != 0) {
        return 1;
    }

    // 主流程：输入校验 -> 预处理 -> 比对 -> 导出
    return runMainPipeline(common_args, argc, argv);
}
