// RaMAx.cpp: 定义应用程序的入口点。
// 主程序：负责解析命令行参数，初始化配置，执行基因组比对流程

#include "SeqPro.h"
#include "anchor.h"
#include "config.hpp"
#include "data_process.h"
#include "index.h"
#include "anchor.h"
#include "rare_aligner.h"
#include "sequence_utils.h"
#include <limits>
#include <algorithm>
#include <variant>

// ------------------------------------------------------------------
// 通用命令行参数结构体（可序列化）
// ------------------------------------------------------------------
struct CommonArgs {
    std::filesystem::path reference_path = "";
    std::filesystem::path query_path = "";
    std::filesystem::path output_path = "";
    std::filesystem::path work_dir_path = "";

    uint_t chunk_size = 10000000;
    uint_t overlap_size = 100000;
    uint_t min_anchor_length = 20;
    uint_t max_anchor_frequency = 50;
    bool restart = false;
    int thread_num = std::thread::hardware_concurrency(); // 默认使用所有 CPU 核心
    PairGenomeOutputFormat output_format = PairGenomeOutputFormat::UNKNOWN;
    bool enable_repeat_masking = false; // 是否启用重复序列遮蔽，默认为 false

    // 新增算法参数
    SearchMode search_mode = FAST_SEARCH;      // 搜索模式
    bool allow_MEM = false;                    // 是否允许MEM
    bool fast_build = true;                    // 是否使用快速索引构建
    SeqPro::Length sampling_interval = 32;    // 采样间隔
    uint_t min_span = 50;                      // 最小跨度阈值

    // 图构建算法选择
    std::string graph_algorithm = "greedy";   // 图构建算法

    // 输出控制参数
    bool include_coordinates = true;           // MAF输出包含坐标
    bool include_sequences = true;             // MAF输出包含序列
    bool compress_output = false;              // 压缩输出

    // 日志相关参数
    std::string log_level = "info";            // 日志级别
    bool verbose = false;                      // 详细输出模式
    bool quiet = false;                        // 静默模式

    // 质量控制参数
    float min_identity = 0.0;                  // 最小序列相似度
    uint_t max_gaps = 0;                       // 最大gap数量
    uint_t filter_short = 0;                   // 过滤短比对
    uint_t min_cluster_size = 0;               // 最小cluster大小

    // 验证和调试参数
    bool validate_anchors = false;             // 验证锚点
    bool validate_clusters = false;            // 验证clusters
    bool validate_graph = false;               // 验证图结构
    bool keep_temp = false;                    // 保留临时文件

    // 支持 cereal 序列化
    template<class Archive>
    void serialize(Archive &ar) {
        ar(
            CEREAL_NVP(reference_path),
            CEREAL_NVP(query_path),
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
            CEREAL_NVP(graph_algorithm),
            CEREAL_NVP(include_coordinates),
            CEREAL_NVP(include_sequences),
            CEREAL_NVP(compress_output),
            CEREAL_NVP(log_level),
            CEREAL_NVP(verbose),
            CEREAL_NVP(quiet),
            CEREAL_NVP(min_identity),
            CEREAL_NVP(max_gaps),
            CEREAL_NVP(filter_short),
            CEREAL_NVP(min_cluster_size),
            CEREAL_NVP(validate_anchors),
            CEREAL_NVP(validate_clusters),
            CEREAL_NVP(validate_graph),
            CEREAL_NVP(keep_temp)
        );
    }
};



// ------------------------------------------------------------------
// 美化参数输出函数
// ------------------------------------------------------------------
inline void printRunConfiguration(const CommonArgs &args) {
    spdlog::info("");
    spdlog::info("============================================================");
    spdlog::info("                     RUN CONFIGURATION                     ");
    spdlog::info("============================================================");
    
    // Input/Output section
    spdlog::info("Input/Output:");
    spdlog::info("  Reference file   : {}", args.reference_path.string());
    spdlog::info("  Query file       : {}", args.query_path.string());
    spdlog::info("  Output file      : {}", args.output_path.string());
    spdlog::info("  Work directory   : {}", args.work_dir_path.string());
    spdlog::info("  Output format    : {}", 
                 args.output_format == PairGenomeOutputFormat::SAM ? "SAM" :
                 args.output_format == PairGenomeOutputFormat::MAF ? "MAF" :
                 args.output_format == PairGenomeOutputFormat::PAF ? "PAF" : "Unknown");
    
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
    spdlog::info("  Min span              : {}", args.min_span);
    spdlog::info("  Graph algorithm       : {}", args.graph_algorithm);
    spdlog::info("  Repeat masking        : {}", args.enable_repeat_masking ? "Enabled" : "Disabled");
    
    spdlog::info("");
    
    // Quality control section
    spdlog::info("Quality Control:");
    if (args.min_identity > 0.0) {
        spdlog::info("  Min identity          : {:.2f}", args.min_identity);
    }
    if (args.max_gaps > 0) {
        spdlog::info("  Max gaps              : {}", args.max_gaps);
    }
    if (args.filter_short > 0) {
        spdlog::info("  Filter short          : {}", args.filter_short);
    }
    if (args.min_cluster_size > 0) {
        spdlog::info("  Min cluster size      : {}", args.min_cluster_size);
    }
    
    spdlog::info("");
    
    // Performance section
    spdlog::info("Performance:");
    spdlog::info("  Thread count          : {}", args.thread_num);
    spdlog::info("  Restart mode          : {}", args.restart ? "Enabled" : "Disabled");
    
    spdlog::info("");
    
    // Output control section
    spdlog::info("Output Control:");
    spdlog::info("  Log level             : {}", args.log_level);
    spdlog::info("  Verbose mode          : {}", args.verbose ? "Enabled" : "Disabled");
    spdlog::info("  Quiet mode            : {}", args.quiet ? "Enabled" : "Disabled");
    spdlog::info("  Include coordinates   : {}", args.include_coordinates ? "Enabled" : "Disabled");
    spdlog::info("  Include sequences     : {}", args.include_sequences ? "Enabled" : "Disabled");
    spdlog::info("  Compress output       : {}", args.compress_output ? "Enabled" : "Disabled");
    spdlog::info("  Keep temp files       : {}", args.keep_temp ? "Enabled" : "Disabled");
    
    spdlog::info("");
    
    // Validation section
    spdlog::info("Validation:");
    spdlog::info("  Validate anchors      : {}", args.validate_anchors ? "Enabled" : "Disabled");
    spdlog::info("  Validate clusters     : {}", args.validate_clusters ? "Enabled" : "Disabled");
    spdlog::info("  Validate graph        : {}", args.validate_graph ? "Enabled" : "Disabled");
    
    spdlog::info("============================================================");
    spdlog::info("");
}

// ------------------------------------------------------------------
// CLI11 参数注册：配置常用命令行参数（参考 RaMAx 主程序）
// ------------------------------------------------------------------
inline void setupCommonOptions(CLI::App *cmd, CommonArgs &args) {
    auto fmt = std::make_shared<CustomFormatter>();
    fmt->column_width(50);
    cmd->formatter(fmt);

    cmd->set_version_flag("-v,--version", std::string("RaMAx version ") + VERSION);

    auto *ref_opt = cmd->add_option("-r,--reference", args.reference_path,
                                    "Path to the reference genome file (FASTA format).")
            ->group("Input Files")
            ->type_name("<path>")->transform(trim_whitespace);

    auto *qry_opt = cmd->add_option("-q,--query", args.query_path,
                                    "Path to the query genome file (FASTA format).")
            ->group("Input Files")
            ->type_name("<path>")->transform(trim_whitespace);

    auto *output_opt = cmd->add_option("-o,--output", args.output_path,
                                       "Path to save alignment results.")
            ->group("Output")->type_name("<path>")
            ->transform(trim_whitespace);

    auto *workspace_opt = cmd->add_option("-w,--workdir", args.work_dir_path,
                                          "Path to the working directory for temporary files.")
            ->group("Output")->type_name("<path>")
            ->transform(trim_whitespace);

    auto *chunk_size_opt = cmd->add_option("--chunk_size", args.chunk_size,
                                           "Size of each chunk for parallel processing (default: 10000000).")
            ->default_val(10000000)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(1000000, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *overlap_size_opt = cmd->add_option("--overlap_size", args.overlap_size,
                                             "Size of overlap between chunks (default: 100000).")
            ->default_val(100000)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *min_anchor_length_opt = cmd->add_option(
                "--min_anchor_length", args.min_anchor_length,
                "Minimum anchor length (default: 20).")
            ->default_val(20)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(1, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *max_anchor_frequency_opt = cmd->add_option(
                "--max_anchor_frequency", args.max_anchor_frequency,
                "Maximum anchor frequency filter (default: 50).")
            ->default_val(50)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    // 新增算法参数
    auto *search_mode_opt = cmd->add_option("--search-mode", args.search_mode,
                                            "Anchor search mode: fast/middle/accurate (default: fast).")
            ->default_val(FAST_SEARCH)
            ->capture_default_str()
            ->group("Software Parameters")
            ->type_name("<mode>")
            ->transform(CLI::CheckedTransformer(
                std::map<std::string, SearchMode>{
                    {"fast", FAST_SEARCH},
                    {"middle", MIDDLE_SEARCH},
                    {"accurate", ACCURATE_SEARCH}
                }, CLI::ignore_case));

    auto *allow_mem_flag = cmd->add_flag("--allow-mem", args.allow_MEM,
                                         "Allow MEM (Maximal Exact Match) instead of only MUM.")
            ->group("Software Parameters");

    auto *slow_build_flag = cmd->add_flag("--slow-build",
                                          "Use slow but more accurate index building method.")
            ->group("Software Parameters");

    auto *sampling_interval_opt = cmd->add_option("--sampling-interval", args.sampling_interval,
                                                  "Reference sequence sampling interval (default: 32).")
            ->default_val(32)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(1, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *min_span_opt = cmd->add_option("--min-span", args.min_span,
                                         "Minimum span threshold for graph construction (default: 50).")
            ->default_val(50)
            ->capture_default_str()
            ->group("Software Parameters")
            ->check(CLI::Range(1, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *graph_algorithm_opt = cmd->add_option("--graph-algorithm", args.graph_algorithm,
                                                "Graph construction algorithm: greedy/dp (default: greedy).")
            ->default_val("greedy")
            ->capture_default_str()
            ->group("Software Parameters")
            ->type_name("<algo>")
            ->transform(CLI::CheckedTransformer(
                std::map<std::string, std::string>{
                    {"greedy", "greedy"},
                    {"dp", "dp"}
                }, CLI::ignore_case));

    // Repeat masking is currently disabled because it depends on external
    // windowmasker binaries that are no longer bundled with RaMAx.
    // cmd->add_flag("--mask-repeats", args.enable_repeat_masking,
    //               "Enable repeat sequence masking.")
    //         ->group("Software Parameters");

    // 质量控制参数
    auto *min_identity_opt = cmd->add_option("--min-identity", args.min_identity,
                                             "Minimum sequence identity threshold (0.0-1.0).")
            ->group("Quality Control")
            ->check(CLI::Range(0.0, 1.0))
            ->type_name("<float>")
            ->transform(trim_whitespace);

    auto *max_gaps_opt = cmd->add_option("--max-gaps", args.max_gaps,
                                         "Maximum number of gaps allowed.")
            ->group("Quality Control")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *filter_short_opt = cmd->add_option("--filter-short", args.filter_short,
                                             "Filter alignments shorter than specified length.")
            ->group("Quality Control")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    auto *min_cluster_size_opt = cmd->add_option("--min-cluster-size", args.min_cluster_size,
                                                 "Minimum cluster size threshold.")
            ->group("Quality Control")
            ->check(CLI::Range(0, std::numeric_limits<int>::max()))
            ->type_name("<int>")
            ->transform(trim_whitespace);

    // 性能参数
    auto *threads_opt = cmd->add_option("-t,--threads", args.thread_num,
                                        "Number of threads to use for parallel processing (default: system cores).")
            ->default_val(std::thread::hardware_concurrency())
            ->envname("RAMA_G_THREADS")
            ->capture_default_str()
            ->group("Performance")
            ->check(CLI::Range(1, std::numeric_limits<int>::max()))
            ->type_name("<int>")->transform(trim_whitespace);

    auto *restart_flag = cmd->add_flag("--restart", args.restart,
                                       "Restart the alignment process by skipping the existing index files.")
            ->group("Performance");

    // 输出控制参数
    auto *no_coordinates_flag = cmd->add_flag("--no-coordinates",
                                              "Do not include coordinates in MAF output.")
            ->group("Output Control");

    auto *no_sequences_flag = cmd->add_flag("--no-sequences",
                                            "Do not include sequences in MAF output.")
            ->group("Output Control");

    auto *compress_output_flag = cmd->add_flag("--compress-output", args.compress_output,
                                               "Compress output files.")
            ->group("Output Control");

    // 日志和输出控制参数
    auto *log_level_opt = cmd->add_option("--log-level", args.log_level,
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
                }, CLI::ignore_case));

    auto *verbose_flag = cmd->add_flag("--verbose", args.verbose,
                                       "Enable verbose output mode.")
            ->group("Output Control");

    auto *quiet_flag = cmd->add_flag("--quiet", args.quiet,
                                     "Enable quiet mode (only errors).")
            ->group("Output Control");

    auto *keep_temp_flag = cmd->add_flag("--keep-temp", args.keep_temp,
                                         "Keep temporary files for debugging.")
            ->group("Output Control");

    // 验证和调试参数
    auto *validate_anchors_flag = cmd->add_flag("--validate-anchors", args.validate_anchors,
                                                "Validate anchor sequence matching correctness.")
            ->group("Validation");

    auto *validate_clusters_flag = cmd->add_flag("--validate-clusters", args.validate_clusters,
                                                 "Validate cluster correctness.")
            ->group("Validation");

    auto *validate_graph_flag = cmd->add_flag("--validate-graph", args.validate_graph,
                                              "Validate graph structure correctness.")
            ->group("Validation");

    // Set dependencies and exclusions
    restart_flag->needs(workspace_opt);
    restart_flag->excludes(ref_opt,
                           qry_opt, output_opt, threads_opt,
                           chunk_size_opt, overlap_size_opt,
                           min_anchor_length_opt, max_anchor_frequency_opt,
                           search_mode_opt, allow_mem_flag, slow_build_flag,
                           sampling_interval_opt, min_span_opt, graph_algorithm_opt,
                           min_identity_opt, max_gaps_opt, filter_short_opt,
                           min_cluster_size_opt);

    // 互斥选项
    verbose_flag->excludes(quiet_flag);
    quiet_flag->excludes(verbose_flag);

    // 处理特殊标志的逻辑将在main函数中处理
}


int main(int argc, char **argv) {
    // 初始化异步日志线程池（spdlog）
    spdlog::init_thread_pool(8192, 1); // 日志缓冲区容量 8192，单线程日志写入
    // setupLogger();  // 控制台输出日志（不带文件）
#ifdef _DEBUG_
    spdlog::set_level(spdlog::level::debug);  // 让 debug 级别可见
#endif


    // 初始化 CLI 命令行应用
    CLI::App app{"RaMA-G: A High-performance Genome Alignment Tool"};

    CommonArgs common_args; // 存储用户输入的参数
    setupCommonOptions(&app, common_args); // 注册参数解析
    CLI11_PARSE(app, argc, argv); // 开始解析命令行参数

    // 处理特殊标志的逻辑
    if (app.count("--slow-build")) {
        common_args.fast_build = false;
    }
    if (app.count("--no-coordinates")) {
        common_args.include_coordinates = false;
    }
    if (app.count("--no-sequences")) {
        common_args.include_sequences = false;
    }

    // 设置日志级别
    if (common_args.quiet) {
        spdlog::set_level(spdlog::level::err);
    } else if (common_args.verbose) {
        spdlog::set_level(spdlog::level::debug);
    } else {
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

    try {
        // ------------------------------
        // 模式 1：重启模式（--restart）
        // ------------------------------
        if (common_args.restart) {
            if (common_args.work_dir_path.empty()) {
                throw CLI::RequiredError("In restart mode, --workdir (-w) is required.");
            }

            // 检查工作目录是否存在
            if (std::filesystem::exists(common_args.work_dir_path)) {
                if (!std::filesystem::is_directory(common_args.work_dir_path)) {
                    throw CLI::ValidationError("Work directory is not valid: " + common_args.work_dir_path.string());
                }
            } else {
                std::filesystem::create_directories(common_args.work_dir_path);
            }

            // 初始化日志器
            setupLoggerWithFile(common_args.work_dir_path);
            spdlog::info("RaMA-G version {}", VERSION);
#ifdef _DEBUG_
            spdlog::info("Restart mode enabled.");
#endif

            // 加载之前保存的参数配置文件
            FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
            std::ifstream is(config_path);
            if (!is) {
                spdlog::error("Failed to open {} for loading CommonArgs", config_path.string());
                return false;
            }
            cereal::JSONInputArchive archive(is);
            archive(common_args); // 反序列化参数
#ifdef _DEBUG_
            spdlog::info("CommonArgs loaded from {}", config_path.string());
#endif
        }

        // ------------------------------
        // 模式 2：正常运行模式
        // ------------------------------
        else {
            // 检查必要参数
            if (common_args.reference_path.empty())
                throw CLI::RequiredError("Missing required option: --reference (-r)");
            if (common_args.query_path.empty())
                throw CLI::RequiredError("Missing required option: --query (-q)");
            if (common_args.output_path.empty())
                throw CLI::RequiredError("Missing required option: --output (-o)");
            if (common_args.work_dir_path.empty())
                throw CLI::RequiredError("Missing required option: --workdir (-w)");

#ifndef _DEBUG_
			// 非调试模式下：确保工作目录为空
			if (std::filesystem::exists(common_args.work_dir_path)) {
				if (!std::filesystem::is_directory(common_args.work_dir_path)) {
					throw CLI::ValidationError("Work directory is not valid: " + common_args.work_dir_path.string());
				}
				if (!std::filesystem::is_empty(common_args.work_dir_path)) {
					throw CLI::ValidationError("Work directory is not empty: " + common_args.work_dir_path.string());
				}
			}
			else {
				std::filesystem::create_directories(common_args.work_dir_path);
			}
#endif

            setupLoggerWithFile(common_args.work_dir_path);
            spdlog::info("RaMA-G version {}", VERSION);

            // 验证参考文件路径
            std::string ref_str = common_args.reference_path.string();
            if (isUrl(ref_str)) {
                verifyUrlReachable(ref_str);
            } else {
                verifyLocalFile(common_args.reference_path);
            }

            // 验证查询文件路径
            std::string qry_str = common_args.query_path.string();
            if (isUrl(qry_str)) {
                verifyUrlReachable(qry_str);
            } else {
                verifyLocalFile(common_args.query_path);
            }

            // 自动识别输出格式（支持 .sam / .maf / .paf）
            common_args.output_format = detectPairGenomeOutputFormat(common_args.output_path);
            if (common_args.output_format == PairGenomeOutputFormat::UNKNOWN) {
                throw std::runtime_error("Invalid output file extension. Supported: .sam, .maf, .paf");
            }

            // 检查 chunk 和 overlap 参数
            if (common_args.overlap_size >= common_args.chunk_size) {
                throw std::runtime_error("Overlap size must be less than chunk size.");
            }

            // 保存参数配置文件（用于 --restart）
            FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
            std::ofstream os(config_path);
            if (!os) {
                spdlog::error("Failed to open {} for saving CommonArgs", config_path.string());
                return false;
            }
            cereal::JSONOutputArchive archive(os);
            archive(cereal::make_nvp("common_args", common_args));
#ifdef _DEBUG_
            spdlog::info("Configuration saved to {}", config_path.string());
#endif
        }
    } catch (const std::runtime_error &e) {
        spdlog::error("{}", e.what());
        spdlog::error("Use --help for usage information.");
        spdlog::error("Exiting with error code 1.");
        return 1;
    }

    // 显示运行配置
    printRunConfiguration(common_args);

    // ------------------------------
    // 主流程开始
    // ------------------------------
    spdlog::info("Executed command: {}", getCommandLine(argc, argv));
    spdlog::info("");
    spdlog::info("============================================================");
    spdlog::info("                      INPUT VALIDATION                     ");
    spdlog::info("============================================================");
    
    spdlog::info("Reference genome: {} (size: {})",
                 common_args.reference_path.string(),
                 getReadableFileSize(common_args.reference_path));

    spdlog::info("Query genome: {} (size: {})",
                 common_args.query_path.string(),
                 getReadableFileSize(common_args.query_path));

    spdlog::info("");
    spdlog::info("============================================================");
    spdlog::info("                    DATA PREPROCESSING                     ");
    spdlog::info("============================================================");
    
    SpeciesPathMap species_path_map;
    species_path_map["reference"] = common_args.reference_path;
    species_path_map["query"] = common_args.query_path;

    // 拷贝或下载原始文件（并行执行）
    copyRawData(common_args.work_dir_path, species_path_map, common_args.thread_num);

    // 声明interval文件映射和SeqPro managers变量
    std::map<SpeciesName, FilePath> interval_files_map;
    std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers;
    SeqPro::Length reference_min_seq_length = std::numeric_limits<SeqPro::Length>::max();
    // 清洗 FASTA 文件（统一格式，替换非法字符）
    cleanRawDataset(common_args.work_dir_path, species_path_map,
                    common_args.thread_num);
    if (common_args.enable_repeat_masking) {
        spdlog::warn("Repeat masking is currently disabled; continuing without repeat masking.");
        common_args.enable_repeat_masking = false;
    }
    // ------------------------------
    // 重复遮蔽
    // ------------------------------
    if (common_args.enable_repeat_masking) {
        spdlog::info("Repeat masking enabled. Generating interval files based on "
            "raw files...");

        // 1. 生成interval文件
        interval_files_map = repeatSeqMasking(
            common_args.work_dir_path, species_path_map, common_args.thread_num);

        if (interval_files_map.empty()) {
            // 没成功生成，中止程序建议用户检查一下
            spdlog::error("No interval files were generated. Please check the input "
                "FASTA files and ensure they are valid.");
            return 1;
        }
        spdlog::info("Interval files generated successfully.");
        spdlog::info("Creating SeqPro managers with repeat masking support...");

        for (const auto &[species_name, cleaned_fasta_path]: species_path_map) {
            if (interval_files_map.contains(species_name)) {
                try {
                    auto original_manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
                    auto manager = std::make_unique<SeqPro::MaskedSequenceManager>(
                        std::move(original_manager),
                        interval_files_map[species_name]
                    );
                    #ifdef _DEBUG_
                    spdlog::info("[{}] SeqPro Manager created with repeat masking: {}", 
                                 species_name, cleaned_fasta_path.string());
#endif
                    auto mgr_var = std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
                    seqpro_managers[species_name] = mgr_var;
                    SequenceUtils::recordReferenceSequenceStats(
                        species_name, seqpro_managers[species_name], reference_min_seq_length);
                    // 记录序列统计信息
                    // SequenceUtils::recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);
                    
                    // 移动到seqpro_managers中
                    //seqpro_managers[species_name] = std::move(manager);
                } catch (const std::exception &e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                                  e.what());
                    return 1;
                }
            } else {
                try {
                    auto manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);

#ifdef _DEBUG_
                    spdlog::info("[{}] SeqPro Manager created without repeat masking: {}", 
                                 species_name, cleaned_fasta_path.string());
#endif
                    auto mgr_var = std::make_shared<SeqPro::ManagerVariant>(std::move(manager));
                    seqpro_managers[species_name] = mgr_var;
                    SequenceUtils::recordReferenceSequenceStats(
                        species_name, seqpro_managers[species_name], reference_min_seq_length);
                    // 记录序列统计信息
                    // SequenceUtils::recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);

                     //seqpro_managers[species_name] = std::move(manager);
                } catch (const std::exception &e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                                  e.what());
                    return 1;
                }
            }
        }
    } else {
        // 不开启重复遮蔽，基于清洗后的文件创建常规SeqPro managers
        // 不开启重复遮蔽，基于清洗后的文件创建常规 SeqPro managers
        spdlog::info("Repeat masking disabled. Creating standard SeqPro managers...");

        for (const auto& [species_name, cleaned_fasta_path] : species_path_map) {
            try {
                /* 1. 先建普通 SequenceManager（unique_ptr） */
                auto raw_ptr = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
#ifdef _DEBUG_
                spdlog::info("[{}] SeqPro Manager created: {}", species_name,
                    cleaned_fasta_path.string());
#endif

                /* 2. 把 unique_ptr 包进 ManagerVariant → SharedManagerVariant */
                auto mgr_var = std::make_shared<SeqPro::ManagerVariant>(std::move(raw_ptr));

                /* 3. 存入全局 map */
                seqpro_managers[species_name] = mgr_var;

                /* 4. 记录序列统计（函数需要 SharedManagerVariant） */
                SequenceUtils::recordReferenceSequenceStats(
                    species_name, seqpro_managers[species_name], reference_min_seq_length);

            }
            catch (const std::exception& e) {
                spdlog::error("[{}] Error creating SeqPro manager: {}", species_name, e.what());
                return 1;
            }
        }

    }


    // 可选：按染色体拆分、按 chunk 切割（目前注释掉）
    // splitRawDataToChr(...)
    // splitChrToChunk(...)

    // ------------------------------
    // 初始化比对器
    // ------------------------------
    PairRareAligner pra(
        common_args.work_dir_path,
        common_args.thread_num,
        common_args.chunk_size,
        common_args.overlap_size,
        common_args.min_anchor_length,
        common_args.max_anchor_frequency
    );

    try {
        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                     INDEX BUILDING                        ");
        spdlog::info("============================================================");
        
        auto t_start_build = std::chrono::steady_clock::now();

        pra.buildIndex("reference", *seqpro_managers["reference"], common_args.fast_build); // 可切换 CaPS / divsufsort
        auto t_end_build = std::chrono::steady_clock::now();
        std::chrono::duration<double> build_time = t_end_build - t_start_build;
        spdlog::info("Index built in {:.3f} seconds.", build_time.count());

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                     QUERY ALIGNMENT                       ");
        spdlog::info("============================================================");
        
        auto t_start_align = std::chrono::steady_clock::now();
        sdsl::int_vector<0> ref_global_cache;
        // 初始化ref_global_cache：采样策略避免二分搜索
        auto sampling_interval = std::min(static_cast<SeqPro::Length>(common_args.sampling_interval), reference_min_seq_length);
    
        // 使用工具函数构建缓存
        SequenceUtils::buildRefGlobalCache(seqpro_managers["reference"], sampling_interval, ref_global_cache);

        MatchVec3DPtr anchors = pra.alignPairGenome(
            "query", *seqpro_managers["query"], common_args.search_mode, common_args.allow_MEM,true, ref_global_cache, sampling_interval);
        auto t_end_align = std::chrono::steady_clock::now();
        std::chrono::duration<double> align_time = t_end_align - t_start_align;
        spdlog::info("Query aligned in {:.3f} seconds.", align_time.count());


        // ------------------------------
        // 验证anchors的序列匹配正确性
        // ------------------------------
        if (common_args.validate_anchors) {
            ValidationResult validation_result = validateAnchorsCorrectness(
                anchors, 
                seqpro_managers["reference"], 
                seqpro_managers["query"]
            );
        } 

        RaMesh::RaMeshMultiGenomeGraph graph(seqpro_managers);
        
        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                    ANCHOR FILTERING                       ");
        spdlog::info("============================================================");
        
        auto t_start_filer = std::chrono::steady_clock::now();
        ClusterVecPtrByStrandByQueryRefPtr cluster_vec_ptr = pra.filterPairSpeciesAnchors("query", anchors, *seqpro_managers["query"], graph);
        auto t_end_filer = std::chrono::steady_clock::now();
        std::chrono::duration<double> filter_time = t_end_filer - t_start_filer;
        spdlog::info("Anchors clustered in {:.3f} seconds.", filter_time.count());
        if (common_args.validate_clusters) {
            validateClusters(cluster_vec_ptr);
        }

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                   GRAPH CONSTRUCTION                      ");
        spdlog::info("============================================================");
        
        auto t_start_construct = std::chrono::steady_clock::now();
        // 使用选择的算法构建图
        if (common_args.graph_algorithm == "greedy") {
            pra.constructGraphByGreedy("query", *seqpro_managers["query"], cluster_vec_ptr, graph, common_args.min_span);
        } else if (common_args.graph_algorithm == "dp") {
            // 如果有DP算法的话
            pra.constructGraphByGreedy("query", *seqpro_managers["query"], cluster_vec_ptr, graph, common_args.min_span);
        }

        auto t_end_construct = std::chrono::steady_clock::now();
        std::chrono::duration<double> construct_time = t_end_construct - t_start_construct;
        spdlog::info("Graph constructed in {:.3f} seconds.", construct_time.count());
    
        if (common_args.validate_graph) {
            graph.verifyGraphCorrectness(true);
        }
    
        graph.exportToMaf(common_args.output_path, seqpro_managers, common_args.include_coordinates, common_args.include_sequences);
        
        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                       COMPLETION                          ");
        spdlog::info("============================================================");
    }
    catch (std::exception& e) {
        spdlog::error("Export to result failed: {}", e.what());
        return 1;
    }

    spdlog::info("RaMA-G execution completed successfully!");
    return 0;
}
