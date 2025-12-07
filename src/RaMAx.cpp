// RaMAx.cpp: 定义应用程序的入口点。
// 主程序：负责解析命令行参数，初始化配置，执行多基因组比对流程

#include "SeqPro.h"
#include "data_process.h"
#include "config.hpp"
#include "index.h"
#include "rare_aligner.h"
#include "sequence_utils.h"

// ------------------------------------------------------------------
// 通用命令行参数结构体（可序列化）
// ------------------------------------------------------------------
struct CommonArgs {
    std::filesystem::path input_path = "";
    std::filesystem::path output_path = "";
    std::filesystem::path work_dir_path = "";

    uint_t chunk_size = 10000000;
    uint_t overlap_size = 100000;
    uint_t min_anchor_length = 20;
    uint_t max_anchor_frequency = 50;
    bool restart = false;
    int thread_num = std::thread::hardware_concurrency(); // 默认使用所有 CPU 核心
    MultipleGenomeOutputFormat output_format = MultipleGenomeOutputFormat::UNKNOWN;
    bool enable_repeat_masking = false; // 是否启用重复序列遮蔽，默认为 false

    // 新增算法参数
    SearchMode search_mode = ACCURATE_SEARCH;  // 搜索模式
    bool allow_MEM = false;                    // 是否允许MEM
    bool fast_build = true;                    // 是否使用快速索引构建
    SeqPro::Length sampling_interval = 32;    // 采样间隔
    uint_t min_span = 65;                      // 最小跨度阈值

    // 日志相关参数
    std::string log_level = "info";            // 日志级别
    bool verbose = false;                      // 详细输出模式
    bool quiet = false;                        // 静默模式


    // HAL root name
    std::string root_name = "root";

    // 支持 cereal 序列化
    template<class Archive>
    void serialize(Archive& ar) {
        ar(
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
            CEREAL_NVP(root_name)
        );
    }
};

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
// CLI11 参数注册：配置常用命令行参数（参考 RaMAx 主程序）
// ------------------------------------------------------------------
inline void setupCommonOptions(CLI::App* cmd, CommonArgs& args) {
    auto fmt = std::make_shared<CustomFormatter>();
    fmt->column_width(50);
    cmd->formatter(fmt);

    cmd->set_version_flag("-v,--version", std::string("RaMAx version ") + VERSION);

    auto* input_opt = cmd->add_option("-i,--input", args.input_path,
        "Path to the sequence file (txt format).")
        ->group("Input Files")
        ->type_name("<path>")->transform(trim_whitespace);


    auto* output_opt = cmd->add_option("-o,--output", args.output_path,
        "Path to save alignment results.")
        ->group("Output")->type_name("<path>")
        ->transform(trim_whitespace);

    auto* workspace_opt = cmd->add_option("-w,--workdir", args.work_dir_path,
        "Path to the working directory for temporary files.")
        ->group("Output")->type_name("<path>")
        ->transform(trim_whitespace);

    auto* chunk_size_opt = cmd->add_option("--chunk_size", args.chunk_size,
        "Size of each chunk for parallel processing (default: 10000000).")
        ->default_val(10000000)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1000000, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    auto* root_opt = cmd->add_option("--root", args.root_name,
        "Root genome name used in HAL (default: 'root')")
        ->group("Output")
        ->type_name("<string>")
        ->transform(trim_whitespace);

    auto* overlap_size_opt = cmd->add_option("--overlap_size", args.overlap_size,
        "Size of overlap between chunks (default: 100000).")
        ->default_val(0)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    auto* min_anchor_length_opt = cmd->add_option(
        "--min_anchor_length", args.min_anchor_length,
        "Minimum anchor length (default: 20).")
        ->default_val(20)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    auto* max_anchor_frequency_opt = cmd->add_option(
        "--max_anchor_frequency", args.max_anchor_frequency,
        "Maximum anchor frequency filter (default: 50).")
        ->default_val(50)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(0, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    // 新增算法参数
    auto* search_mode_opt = cmd->add_option("--search-mode", args.search_mode,
        "Anchor search mode: fast/middle/accurate (default: accurate).")
        ->default_val(ACCURATE_SEARCH)
        ->capture_default_str()
        ->group("Software Parameters")
        ->type_name("<mode>")
        ->transform(CLI::CheckedTransformer(
            std::map<std::string, SearchMode>{
                {"fast", FAST_SEARCH},
                { "middle", MIDDLE_SEARCH },
                { "accurate", ACCURATE_SEARCH }
    }, CLI::ignore_case));

    auto* allow_mem_flag = cmd->add_flag("--allow-mem", args.allow_MEM,
        "Allow MEM (Maximal Exact Match) instead of only MUM.")
        ->group("Software Parameters");

    auto* slow_build_flag = cmd->add_flag("--slow-build",
        "Use slow but more accurate index building method.")
        ->group("Software Parameters");

    auto* sampling_interval_opt = cmd->add_option("--sampling-interval", args.sampling_interval,
        "Reference sequence sampling interval (default: 32).")
        ->default_val(32)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    auto* min_span_opt = cmd->add_option("--min-span", args.min_span,
        "Minimum span threshold for graph construction (default: 50).")
        ->default_val(65)
        ->capture_default_str()
        ->group("Software Parameters")
        ->check(CLI::Range(1, std::numeric_limits<int>::max()))
        ->type_name("<int>")
        ->transform(trim_whitespace);

    /* auto* mask_repeats_flag = */
    cmd->add_flag("--mask-repeats", args.enable_repeat_masking,
        "Enable repeat sequence masking.")
        ->group("Software Parameters");


    // 性能参数
    auto* threads_opt = cmd->add_option("-t,--threads", args.thread_num,
        "Number of threads to use for parallel processing (default: system cores).")
        ->default_val(std::thread::hardware_concurrency())
        ->envname("RAMAx_THREADS")
        ->capture_default_str()
        ->group("Performance")
        ->check(CLI::Range(1, std::numeric_limits<int>::max()))
        ->type_name("<int>")->transform(trim_whitespace);



    auto* restart_flag = cmd->add_flag("--restart", args.restart,
        "Restart the alignment process by skipping the existing index files.")
        ->group("Performance");

    // 日志和输出控制参数
    auto* log_level_opt = cmd->add_option("--log-level", args.log_level,
        "Log level: debug/info/warn/error (default: info).")
        ->default_val("info")
        ->capture_default_str()
        ->group("Output Control")
        ->type_name("<level>")
        ->transform(CLI::CheckedTransformer(
            std::map<std::string, std::string>{
                {"debug", "debug"},
                { "info", "info" },
                { "warn", "warn" },
                { "error", "error" }
    }, CLI::ignore_case));

    auto* verbose_flag = cmd->add_flag("--verbose", args.verbose,
        "Enable verbose output mode.")
        ->group("Output Control");

    auto* quiet_flag = cmd->add_flag("--quiet", args.quiet,
        "Enable quiet mode (only errors).")
        ->group("Output Control");



    // Set dependencies and exclusions
    restart_flag->needs(workspace_opt);
    restart_flag->excludes(input_opt,
        output_opt, threads_opt,
        chunk_size_opt, overlap_size_opt,
        min_anchor_length_opt, max_anchor_frequency_opt,
        search_mode_opt, allow_mem_flag, slow_build_flag,
        sampling_interval_opt, min_span_opt, root_opt);

    // 互斥选项
    verbose_flag->excludes(quiet_flag);
    quiet_flag->excludes(verbose_flag);

    // 处理 --slow-build 标志的逻辑
    // 在解析完成后检查slow_build标志并相应设置fast_build
}


int main(int argc, char** argv) {
    // 初始化异步日志线程池（spdlog）
    spdlog::init_thread_pool(8192, 1); // 日志缓冲区容量 8192，单线程日志写入
    // setupLogger();  // 控制台输出日志（不带文件）
    // 初始化 CLI 命令行应用
    CLI::App app{ "RaMAx: A High-performance Genome Alignment Tool" };

    CommonArgs common_args; // 存储用户输入的参数
    setupCommonOptions(&app, common_args); // 注册参数解析
    CLI11_PARSE(app, argc, argv); // 开始解析命令行参数

    // 处理特殊标志的逻辑
    if (app.count("--slow-build")) {
        common_args.fast_build = false;
    }

    // 设置日志级别
    if (common_args.quiet) {
        spdlog::set_level(spdlog::level::err);
    }
    else if (common_args.verbose) {
        spdlog::set_level(spdlog::level::debug);
    }
    else {
        if (common_args.log_level == "debug") {
            spdlog::set_level(spdlog::level::debug);
        }
        else if (common_args.log_level == "info") {
            spdlog::set_level(spdlog::level::info);
        }
        else if (common_args.log_level == "warn") {
            spdlog::set_level(spdlog::level::warn);
        }
        else if (common_args.log_level == "error") {
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
            }
            else {
                std::filesystem::create_directories(common_args.work_dir_path);
            }

            // 初始化日志器
            setupLoggerWithFile(common_args.work_dir_path);
            spdlog::info("RaMAx version {}", VERSION);
            spdlog::info("Restart mode enabled.");

            // 加载之前保存的参数配置文件
            FilePath config_path = common_args.work_dir_path / CONFIG_FILE;
            std::ifstream is(config_path);
            if (!is) {
                spdlog::error("Failed to open {} for loading CommonArgs", config_path.string());
                return false;
            }
            cereal::JSONInputArchive archive(is);
            archive(common_args); // 反序列化参数
            spdlog::info("CommonArgs loaded from {}", config_path.string());
        }

        // ------------------------------
        // 模式 2：正常运行模式
        // ------------------------------
        else {
            // 检查必要参数
            if (common_args.input_path.empty())
                throw CLI::RequiredError("Missing required option: --input (-i)");
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
            spdlog::info("RaMAx version {}", VERSION);
            spdlog::info("Multiple genome alignment mode enabled.");


            common_args.output_format = detectMultipleGenomeOutputFormat(common_args.output_path);
            if (common_args.output_format == MultipleGenomeOutputFormat::UNKNOWN) {
                throw std::runtime_error("Invalid output file extension. Supported: .hal, .maf");
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
            spdlog::info("Configuration saved to {}", config_path.string());
        }
    }
    catch (const std::runtime_error& e) {
        spdlog::error("{}", e.what());
        spdlog::error("Use --help for usage information.");
        spdlog::error("Exiting with error code 1.");
        return 1;
    }

    // 显示运行配置
    printRunConfiguration(common_args);

    //	// ------------------------------
    //	// 主流程开始
    //	// ------------------------------
    spdlog::info("Executed command: {}", getCommandLine(argc, argv));
    spdlog::info("");
    spdlog::info("============================================================");
    spdlog::info("                      INPUT VALIDATION                     ");
    spdlog::info("============================================================");

    SpeciesPathMap species_path_map;
    NewickParser newick_tree;

    try {
        std::string root = common_args.root_name;
        parseSeqfile(
            common_args.input_path,
            newick_tree, species_path_map, root);

        // 验证species_path_map里的每个物种的路径是否合法
        for (const auto& [species, path] : species_path_map) {
            if (isUrl(path.string())) {
                verifyUrlReachable(path.string());
            }
            else {
                verifyLocalFile(path);
            }
            spdlog::info("Input genome: {} (size: {})",
                species,
                getReadableFileSize(path));
        }

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                    DATA PREPROCESSING                     ");
        spdlog::info("============================================================");

        // ------------------------------
        // 数据预处理阶段
        // ------------------------------

        // 拷贝或下载原始文件（并行执行）
        copyRawData(common_args.work_dir_path, species_path_map, common_args.thread_num);

        // 声明interval文件映射和SeqPro managers变量
        std::map<SpeciesName, FilePath> interval_files_map;
        std::map<SpeciesName, SeqPro::SharedManagerVariant> seqpro_managers;
        SeqPro::Length reference_min_seq_length = std::numeric_limits<SeqPro::Length>::max();

        // 清洗 FASTA 文件（统一格式，替换非法字符）
        cleanRawDataset(common_args.work_dir_path, species_path_map, common_args.thread_num);

        // 如果要重复遮蔽
        if (common_args.enable_repeat_masking) {
            spdlog::info("Repeat masking enabled. Generating interval files based on raw files...");

            // 1. 生成interval文件
            interval_files_map = repeatSeqMasking(
                common_args.work_dir_path, species_path_map, common_args.thread_num);

            if (interval_files_map.empty()) {
                // 没成功生成，中止程序建议用户检查一下
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

                        // 移动到seqpro_managers中
                        seqpro_managers[species_name] = std::move(shared_manager);
                    }
                    catch (const std::exception& e) {
                        spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                            e.what());
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
                        spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                            e.what());
                        return 1;
                    }
                }
            }
        }
        else {
            // 不开启重复遮蔽，基于清洗后的文件创建  MaskedSeqPro managers
            spdlog::info("Repeat masking disabled. Creating standard SeqPro managers...");

            for (const auto& [species_name, cleaned_fasta_path] : species_path_map) {
                try {
                    auto original_manager = std::make_unique<SeqPro::SequenceManager>(cleaned_fasta_path);
                    auto manager = std::make_unique<SeqPro::MaskedSequenceManager>(
                        std::move(original_manager)
                    );
                    spdlog::info("[{}] SeqPro Manager created: {}", species_name,
                        cleaned_fasta_path.string());
                    // 记录序列统计信息
                    SequenceUtils::recordReferenceSequenceStats(species_name, manager, reference_min_seq_length);


                    auto shared_manager = std::make_shared<SeqPro::ManagerVariant>(std::move(manager));

                    seqpro_managers[species_name] = std::move(shared_manager);
                }
                catch (const std::exception& e) {
                    spdlog::error("[{}] Error creating SeqPro manager: {}", species_name,
                        e.what());
                    return 1;
                }
            }
        }


        // ------------------------------
        // 初始化比对器
        // ------------------------------
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

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                    STAR ALIGNMENT                         ");
        spdlog::info("============================================================");

        // ------------------------------
        // 步骤 1：构建索引
        // ------------------------------
        auto t_start_align = std::chrono::steady_clock::now();

        // 初始化ref_global_cache：采样策略避免二分搜索
        auto sampling_interval = std::min(static_cast<SeqPro::Length>(common_args.sampling_interval), reference_min_seq_length);
        uint_t tree_root = 0;

        std::unique_ptr<RaMesh::RaMeshMultiGenomeGraph> graph = mra.starAlignment(seqpro_managers, tree_root, common_args.search_mode, common_args.fast_build, common_args.allow_MEM, common_args.enable_repeat_masking, sampling_interval, common_args.min_span);

        auto t_end_align = std::chrono::steady_clock::now();
        std::chrono::duration<double> align_time = t_end_align - t_start_align;

        spdlog::info("");
        spdlog::info("============================================================");
        spdlog::info("                       COMPLETION                          ");
        spdlog::info("============================================================");
        spdlog::info("Star alignment completed in {:.3f} seconds.", align_time.count());

        std::vector<SpeciesName> species_names = newick_tree.getLeafNames();
        // 原始输出文件，比如 "mammals.maf"
        std::filesystem::path out0 = common_args.output_path;

        // 拆分成 parent + stem + ext
        auto parent = out0.parent_path();
        auto stem = out0.stem().string();      // "mammals"
        auto ext = out0.extension().string(); // ".maf"

        // 构建 (SpeciesName, FilePath) 列表
        std::vector<std::pair<SpeciesName, FilePath>> species_maf_files;
        species_maf_files.reserve(species_names.size());

        for (auto const& sp : species_names) {
            // 生成新文件名：  mammals.<sp>.maf
            std::string filename = stem + "." + sp + ext;
            std::filesystem::path p = parent / filename;
            species_maf_files.emplace_back(sp, p);
        }

        // 清理seqpro_managers的所有遮蔽区间
        for (const auto& [species_name, seq_mgr_variant] : seqpro_managers) {
            std::visit([&](const auto& seq_mgr) {
                seq_mgr->clearMaskedRegions();
                }, *seq_mgr_variant);
        }

        // 根据输出格式选择导出方法
        switch (common_args.output_format) {
        case MultipleGenomeOutputFormat::MAF:
            spdlog::info("Exporting to MAF format...");
            // TODO双基因组比对模式后续要改为false，目前只是调试
            // 正常是false，调试用true
            graph->exportToMaf(common_args.output_path, seqpro_managers, true, false);
            /// 导出没有反向链的maf仅供调试使用
            // graph->exportToMafWithoutReverse(common_args.output_path, seqpro_managers, true, false);
            /// 导出多个参考maf
            //graph->exportToMultipleMaf(species_maf_files, seqpro_managers, true, false);
            break;

        case MultipleGenomeOutputFormat::HAL:
            spdlog::info("Exporting to HAL format...");
            // 使用Newick树信息导出HAL格式
            {
                // 直接使用已解析并可能裁剪过的 newick_tree，避免重复读取导致 --root 子树失效
                graph->exportToHal(common_args.output_path, seqpro_managers, newick_tree, true, common_args.root_name);
            }
            break;

        default:
            throw std::runtime_error("Unsupported output format for multiple genome alignment");
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