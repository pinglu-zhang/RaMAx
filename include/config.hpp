#ifndef CONFIG_HPP
#define CONFIG_HPP

// ------------------------------------------------------------------
// 引入头文件：功能模块包括线程池、命令行解析、日志系统、序列化库等
// ------------------------------------------------------------------
#include "threadpool.h"                          // 自定义线程池（可能用于并行加速）
#include "CLI11.hpp"                             // CLI11 命令行解析库（本地版本）

#include "spdlog/spdlog.h"                       // spdlog 主头文件
#include "spdlog/sinks/stdout_color_sinks.h"     // 控制台彩色输出 sink
#include "spdlog/sinks/basic_file_sink.h"        // 文件输出 sink
#include "spdlog/async.h"                        // 异步日志支持

#include <cereal/archives/json.hpp>              // JSON 格式序列化
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>            // 二进制格式支持
#include <cereal/types/array.hpp>

#include "sdsl/int_vector.hpp"              // SDSL 索引支持（本地版本）
#include "sdsl/wt_huff.hpp"
#include "sdsl/util.hpp"
#include "sdsl/suffix_arrays.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <memory>
#include <filesystem>
#include <unordered_map>

// ------------------------------------------------------------------
// 通用配置常量
// ------------------------------------------------------------------
#define VERSION "1.0.2"                   // 版本号
#define LOGGER_NAME "logger"              // 默认日志器名称
#define LOGGER_FILE "RaMAx.log"          // 默认日志文件名
#define CONFIG_FILE "config.json"         // 配置文件路径

// 数据文件和工作路径定义
#define DATA_DIR "data"	
#define RAW_DATA_DIR "raw_data"
#define CLEAN_DATA_DIR "clean_data"
#define SPLIT_CHR_DIR "split_chr"
#define CHUNK_DIR "chunk"
#define CHUNK_MAP_FILE "chunk_map.json"

#define INDEX_DIR "index"
#define RESULT_DIR "result"

// ------------------------------------------------------------------
// 类型别名
// ------------------------------------------------------------------
using FilePath = std::filesystem::path;   // 文件路径
using SpeciesName = std::string;          // 物种名称
using ChrName = std::string;              // 染色体名称

using ChrIndex = uint32_t;
// ------------------------------------------------------------------
// 调试与整数精度配置
// ------------------------------------------------------------------
#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef M64
#define M64 0
#endif

// 根据宏定义切换 32 位或 64 位整数
#if M64
typedef int64_t	int_t;
typedef uint64_t uint_t;
#define PRIdN	PRId64
#define U_MAX	UINT64_MAX
#define I_MAX	INT64_MAX
#define I_MIN	INT64_MIN
#else
typedef int32_t int_t;
typedef uint32_t uint_t;
#define PRIdN	PRId32
#define U_MAX	UINT32_MAX
#define I_MAX	INT32_MAX
#define I_MIN	INT32_MIN
#endif


// TODO 这部分序列化内容放在这里不合适，更改到其他位置
// ------------------------------------------------------------------
// cereal 序列化扩展：支持文件路径与 SDSL 结构
// ------------------------------------------------------------------
namespace cereal {

	// std::filesystem::path 的序列化和反序列化
	template<class Archive>
	void save(Archive& ar, const std::filesystem::path& p) {
		ar(p.string());
	}

	template<class Archive>
	void load(Archive& ar, std::filesystem::path& p) {
		std::string tmp;  ar(tmp);
		p = std::filesystem::path(tmp);
	}

	// 将 SDSL 对象序列化为字符串
	template<typename SdslObj>
	inline std::string sdsl_to_string(const SdslObj& obj) {
		std::ostringstream oss(std::ios::binary);
		sdsl::serialize(obj, oss);
		return oss.str();
	}

	// 从字符串反序列化回 SDSL 对象
	template<typename SdslObj>
	inline void string_to_sdsl(const std::string& buf, SdslObj& obj) {
		std::istringstream iss(buf, std::ios::binary);
		sdsl::load(obj, iss);
	}

	// 特化：sdsl::int_vector<0> 的序列化
	template<class Archive>
	void save(Archive& ar, const sdsl::int_vector<0>& v) {
		ar(sdsl_to_string(v));
	}

	template<class Archive>
	void load(Archive& ar, sdsl::int_vector<0>& v) {
		std::string buf;  ar(buf);
		string_to_sdsl(buf, v);
	}

	// 特化：sdsl::wt_huff<bit_vector> 的序列化
	template<class Archive>
	void save(Archive& ar, const sdsl::wt_huff<sdsl::bit_vector>& wt) {
		ar(sdsl_to_string(wt));
	}

	template<class Archive>
	void load(Archive& ar, sdsl::wt_huff<sdsl::bit_vector>& wt) {
		std::string buf;  ar(buf);
		string_to_sdsl(buf, wt);
	}

} // namespace cereal

// ------------------------------------------------------------------
// CLI11 自定义格式器（美化选项输出）
// ------------------------------------------------------------------
class CustomFormatter : public CLI::Formatter {
public:
	CustomFormatter() : Formatter() {}

	// 自定义参数展示样式（带默认值）
	std::string make_option_opts(const CLI::Option* opt) const override {
		if (opt->get_type_size() == 0) return "";
		std::ostringstream out;
		out << " " << opt->get_type_name();
		if (!opt->get_default_str().empty())
			out << " (default: " << opt->get_default_str() << ")";
		return out.str();
	}

	// 提供使用示例
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ramax -i <seqfile.txt> -o <output.{maf|hal}> -w <workdir> [options]\n\n"
			<< "Example:\n"
			<< "  ramax -i seqfile.txt -o output.maf -w workdir -t 8\n\n";
		return out.str();
	}
};

// ------------------------------------------------------------------
// CLI11 自定义 validator：自动去除参数两侧空白
// ------------------------------------------------------------------
inline CLI::Validator trim_whitespace = CLI::Validator(
	[](std::string& s) {
		auto start = s.find_first_not_of(" \t\n\r");
		auto end = s.find_last_not_of(" \t\n\r");
		s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
		return std::string();  // 空字符串表示验证通过
	}, ""
);

// TODO 这是双基因组比对格式的枚举，可能需要根据实际需求进行修改
// ------------------------------------------------------------------
// 枚举：输出格式
// ------------------------------------------------------------------
enum class PairGenomeOutputFormat {
	SAM,
	MAF,
	PAF,
	UNKNOWN
};

// 根据输出文件扩展名自动判断输出格式
inline PairGenomeOutputFormat detectPairGenomeOutputFormat(const std::filesystem::path& output_path) {
	std::string ext = output_path.extension().string();
	if (ext == ".sam") return PairGenomeOutputFormat::SAM;
	if (ext == ".maf") return PairGenomeOutputFormat::MAF;
	if (ext == ".paf") return PairGenomeOutputFormat::PAF;
	return PairGenomeOutputFormat::UNKNOWN;
}

enum class MultipleGenomeOutputFormat {
	HAL,
	MAF,
	UNKNOWN
};

// 根据输出文件扩展名自动判断输出格式
inline MultipleGenomeOutputFormat detectMultipleGenomeOutputFormat(const std::filesystem::path& output_path) {
	std::string ext = output_path.extension().string();
	if (ext == ".hal") return MultipleGenomeOutputFormat::HAL;
	if (ext == ".maf") return MultipleGenomeOutputFormat::MAF;
	return MultipleGenomeOutputFormat::UNKNOWN;
}


// ------------------------------------------------------------------
// 日志系统初始化
// ------------------------------------------------------------------

// 控制台 + 文件日志
inline void setupLoggerWithFile(std::filesystem::path log_dir) {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	std::filesystem::path log_file = log_dir / LOGGER_FILE;
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string(), true);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

	spdlog::sinks_init_list sinks = { console_sink, file_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
	spdlog::set_level(spdlog::level::trace);
	spdlog::flush_every(std::chrono::seconds(3));
}

// 控制台日志
inline void setupLogger() {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	spdlog::sinks_init_list sinks = { console_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
	spdlog::set_level(spdlog::level::trace);
	spdlog::flush_every(std::chrono::seconds(3));
}

// 获取完整命令行字符串（用于日志记录和重现）
inline std::string getCommandLine(int argc, char** argv) {
	std::ostringstream cmd;
	for (int i = 0; i < argc; ++i) {
		cmd << argv[i];
		if (i != argc - 1) cmd << " ";
	}
	return cmd.str();
}

#endif // CONFIG_HPP
