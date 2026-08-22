// RaMAx.cpp: 定义应用程序的入口点。
// 主程序：负责解析命令行参数，初始化配置，执行多基因组比对流程

#include "SeqPro.h"
#include "data_process.h"
#include "config.hpp"
#include "index.h"
#include "minipoa_locator.h"
#include "ramax_version.h"
#include "rare_aligner.h"
#include "sequence_utils.h"
#include "external_msa_runner.h"
#include <cstdlib>

// ------------------------------------------------------------------
// 通用命令行参数结构体（支持 cereal 序列化）
// 用于控制基因组比对 / 搜索 / 组装相关流程的全局参数
// ------------------------------------------------------------------
struct CommonArgs {

    uint32_t schema_version = 1;

    // ========================
    // 输入 / 输出路径相关参数
    // ========================

    std::filesystem::path input_path = "";      // 输入序列文件路径seqfile
    std::filesystem::path output_path = "";     // 最终输出结果路径
    std::filesystem::path work_dir_path = "";   // 工作目录路径（用于中间文件与缓存）


    // ========================
    // 运行控制参数
    // ========================

    bool restart = false;                       // 是否从已有中间结果重新启动任务
    int thread_num = std::thread::hardware_concurrency(); // 使用的线程数（默认使用所有 CPU 核心）

    MultipleGenomeOutputFormat output_format =
        MultipleGenomeOutputFormat::UNKNOWN;    // 多基因组输出格式（如 HAL、MAF 等）

    bool enable_repeat_masking = false;          // 是否启用重复序列遮蔽（Repeat Masking）

    // ========================
    // 搜索 / 算法相关参数
    // ========================

    SearchMode search_mode = ACCURATE_SEARCH;   // 搜索模式（如精确搜索 / 快速搜索）
    bool allow_MEM = false;                     // 是否允许使用 MEM（Maximal Exact Match）
    bool fast_build = true;                     // 是否启用快速索引构建模式
    SeqPro::Length sampling_interval = 32;      // 索引采样间隔（影响速度与内存）
    uint_t min_span = 65;                       // 锚点或匹配的最小跨度阈值

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
            CEREAL_NVP(work_dir_path),
            CEREAL_NVP(chunk_size),
            CEREAL_NVP(overlap_size),
            CEREAL_NVP(restart),
            CEREAL_NVP(thread_num),
            CEREAL_NVP(output_format),
            CEREAL_NVP(enable_repeat_masking),
            CEREAL_NVP(search_mode),
            CEREAL_NVP(allow_MEM),
            CEREAL_NVP(fast_build),
            CEREAL_NVP(sampling_interval),
            CEREAL_NVP(min_span),
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
constexpr uint32_t CONFIG_SCHEMA_VERSION = 1;

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
    spdlog::info("  Output file      : {}", args.output_path.string());
    spdlog::info("  Work directory   : {}", args.work_dir_path.string());
    spdlog::info("  Output format    : {}",
        args.output_format == MultipleGenomeOutputFormat::HAL ? "HAL" :
        args.output_format == MultipleGenomeOutputFormat::MAF ? "MAF" : "Unknown");

    spdlog::info("");

    // Algorithm parameters section
    spdlog::info("Algorithm Parameters:");
    spdlog::info("  Chunk size            : {:L}", args.chunk_size);
    spdlog::info("  Overlap size          : {:L}", args.overlap_size);
    spdlog::info("  Min anchor length     : {}", args.min_anchor_length);
    spdlog::info("  Max anchor frequency  : {}", args.max_anchor_frequency);
    spdlog::info("  Search mode           : {}", SearchModeToString(args.search_mode));
    spdlog::info("  Allow MEM             : {}", args.allow_MEM ? "Enabled" : "Disabled");
    spdlog::info("  Fast build            : {}", args.fast_build ? "Enabled" : "Disabled");
    spdlog::info("  Sampling interval     : {}", args.sampling_interval);
    spdlog::info("  Cluster Min span      : {}", args.min_span);
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
        "Path to a Cactus-compatible seqfile.")
        ->group("Input Files")                      // 帮助信息分组
        ->type_name("<path>")                       // 参数类型显示名称
        ->transform(trim_whitespace);               // 去除首尾空白字符

    // 输出结果路径
    auto* output_opt = cmd->add_option("-o,--output", args.output_path,
        "Output alignment path in MAF or HAL format.")
        ->group("Output")
        ->type_name("<path>")
        ->transform(trim_whitespace);

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
    auto* chunk_size_opt = cmd->add_option("--chunk_size", args.chunk_size,
        "Size of each chunk for parallel processing (default: 10000000).")
        ->default_val(10000000)                     // 默认值
        ->capture_default_str()                     // 在 --help 中显示默认值
        ->group("Software Parameters")
        ->check(CLI::Range(1000000,
            std::numeric_limits<int>::max()))       // 合法取值范围
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // HAL 文件中的根基因组名称
    auto* root_opt = cmd->add_option("--root", args.root_name,
        "Root genome name used in HAL (default: 'root')")
        ->group("Output")
        ->type_name("<string>")
        ->transform(trim_whitespace);

    // 参考基因组名称
    auto* ref_opt = cmd->add_option("--ref", args.ref_name,
        "Ref genome name used in alignment")
        ->group("Software Parameters")
        ->type_name("<string>")
        ->transform(trim_whitespace);

    // 相邻切片之间的重叠长度
    auto* overlap_size_opt = cmd->add_option("--overlap_size", args.overlap_size,
        "Size of overlap between chunks (default: 100000).")
        ->default_val(0)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 锚点的最小长度
    auto* min_anchor_length_opt = cmd->add_option(
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
    auto* max_anchor_frequency_opt = cmd->add_option(
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
    auto* search_mode_opt = cmd->add_option("--search-mode", args.search_mode,
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

    // 是否允许使用 MEM（Maximal Exact Match）
    auto* allow_mem_flag = cmd->add_flag("--allow-mem", args.allow_MEM,
        "Allow MEM (Maximal Exact Match) instead of only MUM.")
        ->group("Software Parameters");

    // 是否仅运行一轮比对流程
    auto* one_round_flag = cmd->add_flag("--one-round", args.one_round,
        "Only run one round for alignment.")
        ->group("Software Parameters");

    auto* optimize_blocks_flag = cmd->add_flag(
        "--optimize-blocks",
        "Explicitly enable the default Block optimizations (enabled by default).")
        ->group("Graph Optimization");

    auto* merge_exact_blocks_flag = cmd->add_flag(
        "--merge-blocks",
        args.merge_exact_contiguous_blocks,
        "Merge compatible neighboring Blocks.")
        ->group("Graph Optimization");

    auto* merge_query_gap_opt = cmd->add_option(
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

    auto* realign_missing_species_span_opt = cmd->add_option(
        "--realign-span",
        args.species_mismatch_realign_max_span,
        "Maximum span for missing-sequence realignment (bp).")
        ->default_val(3000)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->needs(realign_missing_species_flag)
        ->check(CLI::Range(1, 10000));

    auto* zero_gap_merge_span_opt = cmd->add_option(
        "--zero-gap-span",
        args.species_mismatch_zero_gap_max_span,
        "Maximum span for zero-gap missing windows (bp).")
        ->default_val(200)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->check(CLI::Range(1, 3000));

    auto* structural_break_repair_flag = cmd->add_flag(
        "--repair-breaks",
        args.repair_structural_breaks,
        "Repair high-confidence structural discontinuities.")
        ->group("Graph Optimization");

    auto* structural_break_span_opt = cmd->add_option(
        "--break-span",
        args.structural_break_max_span,
        "Maximum structural-break repair span (bp).")
        ->default_val(1000)
        ->capture_default_str()
        ->group("Graph Optimization")
        ->type_name("<bp>")
        ->check(CLI::Range(1, 1000));

    auto* short_block_repair_flag = cmd->add_flag(
        "--merge-short-blocks",
        args.repair_short_blocks,
        "Try to merge short Blocks with banded KSW2.")
        ->group("Graph Optimization");

    // 使用慢但更精确的索引构建方式
    auto* slow_build_flag = cmd->add_flag("--slow-build",
        "Use slow but more accurate index building method.")
        ->group("Software Parameters");

    // 参考序列索引采样间隔
    auto* sampling_interval_opt = cmd->add_option("--sampling-interval",
        args.sampling_interval,
        "Reference sequence sampling interval (default: 32).")
        ->default_val(32)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 构图或链式连接的最小跨度阈值
    auto* min_span_opt = cmd->add_option("--min-span", args.min_span,
        "Minimum span threshold for graph construction (default: 50).")
        ->default_val(65)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // Repeat masking is currently disabled because it depends on external
    // windowmasker binaries that are no longer bundled with RaMAx.
    // cmd->add_flag("--mask-repeats", args.enable_repeat_masking,
    //     "Enable repeat sequence masking.")
    //     ->group("Software Parameters");

    // ========================
    // 性能相关参数
    // ========================

    // 并行线程数
    auto* threads_opt = cmd->add_option("-t,--threads", args.thread_num,
        "Number of threads to use for parallel processing (default: system cores).")
        ->default_val(std::thread::hardware_concurrency())
        ->envname("RAMAx_THREADS")                  // 支持环境变量设置
        ->capture_default_str()
        ->group("Performance")
        ->check(CLI::Range(1,
            std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 是否从已有索引和中间结果重新启动
    auto* restart_flag = cmd->add_flag("--restart", args.restart,
        "Restart the alignment process by skipping the existing index files.")
        ->group("Performance");

    // ========================
    // 日志与输出控制参数
    // ========================

    // 日志级别设置
    auto* log_level_opt = cmd->add_option("--log-level", args.log_level,
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

    // --restart 模式下，禁止使用以下参数
    restart_flag->excludes(
        input_opt,
        output_opt,
        threads_opt,
        chunk_size_opt,
        overlap_size_opt,
        min_anchor_length_opt,
        max_anchor_frequency_opt,
        search_mode_opt,
        allow_mem_flag,
        slow_build_flag,
        sampling_interval_opt,
        min_span_opt,
        root_opt,
        ref_opt,
        optimize_blocks_flag,
        merge_exact_blocks_flag,
        merge_query_gap_opt,
        realign_missing_species_flag,
        realign_missing_species_span_opt,
        zero_gap_merge_span_opt,
        structural_break_repair_flag,
        structural_break_span_opt,
        short_block_repair_flag,
        one_round_flag,
        log_level_opt
    );

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
static int runRestartMode(CommonArgs& common_args) {
    // restart 模式下 workdir 必须提供
    if (common_args.work_dir_path.empty()) {
        throw CLI::RequiredError("In restart mode, --workdir (-w) is required.");
    }

    // 检查/创建工作目录
    ensureWorkDirValidOrCreate(common_args.work_dir_path);

    // 初始化日志器（输出到文件）
    setupLoggerWithFile(common_args.work_dir_path);
    configureLogLevel(common_args);
    spdlog::info("RaMAx version {}", RAMAX_VERSION);
    logBuildMode();
    spdlog::info("Restart mode enabled.");

    // 加载之前保存的参数配置文件
    FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
    std::ifstream is(config_path);
    if (!is) {
        spdlog::error("Failed to open {} for loading CommonArgs", config_path.string());
        return 1;
    }

    try {
        cereal::JSONInputArchive archive(is);
        archive(cereal::make_nvp("common_args", common_args));
    } catch (const std::exception& error) {
        spdlog::error(
            "Restart configuration is incompatible with RaMAx {}: {}",
            RAMAX_VERSION, error.what());
        return 1;
    }
    if (common_args.schema_version != CONFIG_SCHEMA_VERSION) {
        spdlog::error(
            "Unsupported restart schema version {} (expected {}). Old work "
            "directories are not compatible with RaMAx {}.",
            common_args.schema_version, CONFIG_SCHEMA_VERSION,
            RAMAX_VERSION);
        return 1;
    }
    common_args.restart = true;
    configureLogLevel(common_args);
    spdlog::info("CommonArgs loaded from {}", config_path.string());
    return 0;
}

// ------------------------------
// 模式 2：正常运行模式
// - 检查必需参数 input/output/workdir
// - 非 debug 下：工作目录必须为空（否则报错）
// - 保存参数到 workdir/CONFIG_FILE（用于 --restart）
// ------------------------------
static int runNormalMode(CommonArgs& common_args) {
    // 检查必要参数
    if (common_args.input_path.empty())
        throw CLI::RequiredError("Missing required option: --input (-i)");
    if (common_args.output_path.empty())
        throw CLI::RequiredError("Missing required option: --output (-o)");
    if (common_args.work_dir_path.empty())
        throw CLI::RequiredError("Missing required option: --workdir (-w)");

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

    // 根据输出文件扩展名检测输出格式（.hal / .maf）
    common_args.output_format = detectMultipleGenomeOutputFormat(common_args.output_path);
    if (common_args.output_format == MultipleGenomeOutputFormat::UNKNOWN) {
        throw std::runtime_error("Invalid output file extension. Supported: .hal, .maf");
    }

    // 检查 chunk 与 overlap 的合法性
    if (common_args.overlap_size >= common_args.chunk_size) {
        throw std::runtime_error("Overlap size must be less than chunk size.");
    }

    common_args.schema_version = CONFIG_SCHEMA_VERSION;
    FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
    std::ofstream os(config_path);
    if (!os) {
        spdlog::error("Failed to open {} for saving CommonArgs", config_path.string());
        return 1;
    }

    cereal::JSONOutputArchive archive(os);
    archive(cereal::make_nvp("common_args", common_args));
    spdlog::info("Configuration saved to {}", config_path.string());

    return 0;
}

// ------------------------------
// 运行前的配置阶段：
// - 根据 restart 标志选择模式
// - 捕获 runtime_error 并统一打印提示
// ------------------------------
static int prepareRun(CommonArgs& common_args) {
    try {
        // 模式 1：重启模式（--restart）
        if (common_args.restart) {
            return runRestartMode(common_args);
        }

        // 模式 2：正常运行模式
        return runNormalMode(common_args);
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
    const CommonArgs& common_args,
    NewickParser& newick_tree,
    SpeciesPathMap& species_path_map
) {
    // 解析 seqfile：构建 Newick 树与物种-路径映射
    std::string root = common_args.root_name;
    parseSeqfile(common_args.input_path, newick_tree, species_path_map, root);

    // 逐个校验输入基因组路径是否合法（URL 可达 / 本地文件存在）
    for (const auto& [species, path] : species_path_map) {
        if (isUrl(path.string())) {
            verifyUrlReachable(path.string());
        } else {
            verifyLocalFile(path);
        }

        spdlog::info("Input genome: {} (size: {})",
                     species,
                     getReadableFileSize(path));
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
    copyRawData(common_args.work_dir_path, species_path_map, common_args.thread_num);

    // interval 文件映射（重复遮蔽用）
    std::map<SpeciesName, FilePath> interval_files_map;

    // 参考序列最短长度（用于后续 sampling_interval 截断）
    reference_min_seq_length = std::numeric_limits<SeqPro::Length>::max();

    const bool align_softmasked_regions =
        std::getenv("RAMAX_ALIGN_SOFTMASK") != nullptr;
    // HAL always records the original lowercase runs for lossless export.
    // The experimental alignment switch additionally loads those runs as
    // seed-mask intervals while retaining the complete original coordinate
    // space for extension and export.
    if (common_args.output_format == MultipleGenomeOutputFormat::HAL ||
        align_softmasked_regions) {
        cleanRawDatasetWithSoftMaskIndex(
            common_args.work_dir_path, species_path_map,
            softmask_path_map, common_args.thread_num);
    } else {
        cleanRawDataset(
            common_args.work_dir_path, species_path_map,
            common_args.thread_num);
    }

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
    NewickParser& newick_tree,
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
        newick_tree,
        common_args.thread_num,
        common_args.chunk_size,
        common_args.overlap_size,
        common_args.min_anchor_length,
        common_args.max_anchor_frequency
    );
    mra.allow_mem = common_args.allow_MEM;

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
// 导出结果（MAF / HAL）
// ------------------------------
static void exportResults(
    const CommonArgs &common_args, const NewickParser &newick_tree,
    std::map<SpeciesName, SeqPro::SharedManagerVariant> &seqpro_managers,
    const SoftMask::PathMap &softmask_path_map,
    RaMesh::RaMeshMultiGenomeGraph *graph) {
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

  // 根据输出格式选择导出方法
  switch (common_args.output_format) {
  case MultipleGenomeOutputFormat::MAF:
    spdlog::info("Exporting to MAF format...");
    graph->exportToMaf(common_args.output_path, seqpro_managers, false);
    break;

  case MultipleGenomeOutputFormat::HAL:
    spdlog::info("Exporting to HAL format...");
    // 使用已解析并可能裁剪过的 newick_tree，避免重复读取导致 --root 子树失效
    graph->exportToHal(common_args.output_path, seqpro_managers, newick_tree,
                       common_args.root_name,
                       static_cast<int>(common_args.thread_num),
                       softmask_path_map);
    break;

  default:
    throw std::runtime_error(
        "Unsupported output format for multiple genome alignment");
  }
  logCrossAnchorInsertionRepairStats();
  RaMesh::Alignment::ExternalMsaRunner::instance().logSummary();
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

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                    STAR ALIGNMENT                         ");
        spdlog::info("============================================================");

        // ------------------------------
        // Star alignment 阶段
        // ------------------------------
        double align_seconds = 0.0;
        std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> graph =
            runStarAlignment(common_args, species_path_map, newick_tree,
                             seqpro_managers, reference_min_seq_length, align_seconds);

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                       COMPLETION                          ");
        spdlog::info("============================================================");
        spdlog::info("Star alignment completed in {:.3f} seconds.", align_seconds);

        // ------------------------------
        // 导出阶段
        // ------------------------------
        exportResults(common_args, newick_tree, seqpro_managers,
                      softmask_path_map, graph.get());

        // ------------------------------
        // 清理工作目录
        // ------------------------------
        {
            std::error_code ec;
            std::filesystem::remove_all(common_args.work_dir_path, ec);
            if (ec) {
                spdlog::warn("Failed to remove work directory: {}", ec.message());
            } else {
                spdlog::info("Work directory removed: {}", common_args.work_dir_path.string());
            }
        }

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

    // Apply output policy before emitting deprecation warnings.
    configureLogLevel(common_args);
    try {
        applyGraphOptimizationOptions(app, common_args);
    } catch (const CLI::Error& error) {
        return app.exit(error);
    }

    // 运行前准备：根据 restart 与否进行目录/参数/配置文件处理
    if (prepareRun(common_args) != 0) {
        return 1;
    }

    // 主流程：输入校验 -> 预处理 -> 比对 -> 导出
    return runMainPipeline(common_args, argc, argv);
}
