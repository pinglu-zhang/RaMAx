#ifndef DATA_PROCESS_H
#define DATA_PROCESS_H

// -----------------------------
// 包含依赖头文件
// -----------------------------
#include <regex>
#include <string>
#include <stdexcept>
#include <cctype>
#include <curl/curl.h>      // 用于 URL 请求
#include <zlib.h>           // 用于处理 gzip 压缩格式
#include "config.hpp"
#include "anchor.h"
#include "kseq.h"           // 用于解析 FASTA 格式
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <cstdint>
#include "softmask_index.h"

class FastaManager;
// 初始化 kseq 使用 gzFile 类型
KSEQ_INIT(gzFile, gzread)

// 物种名到文件路径的映射
using SpeciesPathMap = std::unordered_map<SpeciesName, FilePath>;
using SpeciesFastaManagerMap = std::unordered_map<SpeciesName, FastaManager>;
using ChrIndexMap = std::unordered_map<ChrName, std::size_t>;

// -----------------------------
// GzFileWrapper：RAII 封装 gzopen/gzclose
// -----------------------------
class GzFileWrapper {
public:
    explicit GzFileWrapper(const std::string& filename, const std::string& mode = "r") {
        fp_ = gzopen(filename.c_str(), mode.c_str());
        if (!fp_) {
            throw std::runtime_error("Failed to open FASTA file: " + filename);
        }
    }
    ~GzFileWrapper() {
        if (fp_) {
            gzclose(fp_);
        }
    }

    gzFile get() const { return fp_; }

    // 禁止复制操作（只允许移动或指针管理）
    GzFileWrapper(const GzFileWrapper&) = delete;
    GzFileWrapper& operator=(const GzFileWrapper&) = delete;

private:
    gzFile fp_{ nullptr };
};

// -----------------------------
// kseq_t 的智能指针删除器
// -----------------------------
struct KseqDeleter {
    void operator()(kseq_t* seq) const {
        if (seq) kseq_destroy(seq);
    }
};

using KseqPtr = std::unique_ptr<kseq_t, KseqDeleter>;

// -----------------------------
// FastaProcessor 类：用于读取、清洗、写入 fasta 序列数据
// -----------------------------
class FastaProcessor {
public:
    FilePath fasta_path_;  // FASTA 文件路径
    bool has_n_in_fasta = false;

    // 构造函数与析构函数：采用 RAII，不再需要在析构函数中手动释放资源
    FastaProcessor() = default;
    explicit FastaProcessor(const FilePath& fasta_path)
        : fasta_path_(fasta_path)
    {
        fasta_open(); // 打开 FASTA 文件并初始化解析器
    }
    ~FastaProcessor() = default; // unique_ptr 自动释放资源

    FastaProcessor(FastaProcessor&&) = default;
    FastaProcessor& operator=(FastaProcessor&&) = default;

    // 获取下一个序列记录（返回 false 表示读完）
    bool nextRecord(std::string& header, std::string& sequence);

    // 重置读取状态，从头开始读取
    void reset();

    // 将多个记录拼接为一个序列（用于压缩或合并）
    std::string concatRecords(char terminator = '\0',
        size_t limit = SIZE_MAX);

    // 获取序列文件的基本统计信息
    struct Stats {
        size_t record_count = 0;
        size_t total_bases = 0;
        size_t min_len = SIZE_MAX;
        size_t max_len = 0;
        size_t average_len = 0;
    };
    Stats getStats();

    // 将清洗后的序列写入新的 fasta 文件
    FilePath writeCleanedFasta(const FilePath& output_file, uint64_t line_width = 60);

protected:
    // 内部成员变量：gzip 文件包装和 kseq_t 管理器
    std::unique_ptr<GzFileWrapper> gz_file_wrapper_;
    KseqPtr seq_;

    bool clean_data_{ false }; // 是否启用清洗

    // 清洗序列（全部转大写 + 过滤非法字符）
    void cleanSequence(std::string& seq);

    // 打开 fasta 文件并初始化 kseq_t
    void fasta_open();
};

// -----------------------------
// FastaManager 类：用于索引和管理 fasta 序列数据
// -----------------------------
class FastaManager : public FastaProcessor {
public:
    FilePath fai_path_;    // 可选 .fai 索引文件路径
    ChrIndexMap idx_map;

    // 构造函数
    FastaManager() = default;

    FastaManager(const FilePath& fasta_path, const FilePath& fai_path = FilePath())
        : FastaProcessor(fasta_path), fai_path_(fai_path)
    {
        if (!fai_path.empty()) {
            loadFaiRecords(fai_path);
        }
        idx_map.reserve(fai_records.size());          // 预留桶，避免重复 rehash

        for (std::size_t i = 0; i < fai_records.size(); ++i)
            idx_map.emplace(fai_records[i].seq_name, i);  // 若有重复名称，后插入会被忽略
    }
    FastaManager(const FastaManager&) = delete;
    FastaManager& operator=(const FastaManager&) = delete;

    // 允许移动（unique_ptr 本身可移动，所以 default 就行）
    FastaManager(FastaManager&&)            noexcept = default;
    FastaManager& operator=(FastaManager&&) noexcept = default;

    // 清洗并生成 .fai 索引，返回新 fasta 文件路径
    FilePath cleanAndIndexFasta(const FilePath& output_dir,
        const std::string& prefix,
        uint64_t line_width = 60);

    // -----------------------------
    // FAI 索引记录结构
    // -----------------------------
    struct FaiRecord {
        ChrName seq_name;        // 序列名（如 chr1）
        uint_t global_start_pos;    // 该序列在拼接序列中的全局起始坐标
        uint_t length;              // 碱基总数
        uint_t offset;              // 在文件中的偏移量（跳过注释）
        uint_t line_bases;          // 每行的碱基数
        uint_t line_bytes;          // 每行占用的字节数（含换行）
    };

    // 所有 FAI 记录
    std::vector<FaiRecord> fai_records;

    // 从指定序列中获取区间子串（基于名称 + 坐标）
    std::string getSubSequence(const std::string& seq_name, size_t start, size_t length);

    // 从拼接序列的全局坐标中提取连续片段
    std::string getSubConcatSequence(size_t start, size_t length);

    // 获取拼接后的总序列长度
    uint_t getConcatSeqLength();

    // 从全局坐标范围反查染色体名称
    ChrName getChrName(uint_t global_start, uint_t length);

    // 按 chunk_size 和 overlap_size 对所有序列预分段
    //RegionVec preAllocateChunks(uint_t chunk_size, uint_t overlap_size);

    // 通用的隐藏区间函数，根据提供的区间数据隐藏指定区间，生成新的FASTA文件和对应的FAI索引
    bool hideIntervalsAndGenerateFai(
        const std::map<std::string, std::vector<std::pair<size_t, size_t>>>& intervals_by_seq,
        const FilePath& output_fasta_path,
        const FilePath& output_fai_path,
        size_t line_width = 80);

    // 根据interval文件隐藏指定区间，生成新的FASTA文件和对应的FAI索引
    bool hideIntervalsFromFileAndGenerateFai(
        const FilePath& interval_file_path,
        const FilePath& output_fasta_path,
        const FilePath& output_fai_path,
        size_t line_width = 80);

private:
    // 若不存在 .fai 索引则重新扫描并生成
    bool reScanAndWriteFai(const FilePath& fa_path,
        const FilePath& fai_path,
        size_t line_width) const;

    // 加载已有 .fai 文件
    void loadFaiRecords(const FilePath& fai_path);
};


/**
 * @brief Structure to store information of each node after parsing.
 *
 * @details
 * - id:            Unique identifier for the node (assigned in the order of appearance).
 * - name:          Name of the node (e.g., leaf name or internal node name).
 * - branchLength:  Length of the branch connecting this node to its parent.
 * - isLeaf:        Flag indicating whether the node is a leaf node.
 * - father:        Index of the parent node in the nodes vector, -1 indicates no parent (root node).
 * - leftChild:     Index of the left child node in the nodes vector, -1 indicates no left child.
 * - rightChild:    Index of the right child node in the nodes vector, -1 indicates no right child.
 */
struct NewickTreeNode {
    int      id;
    std::string   name;
    double   branchLength;
    bool     isLeaf;

    int      father;
    int      leftChild;
    int      rightChild;

    NewickTreeNode()
        : id(-1)
        , branchLength(0.0)
        , isLeaf(false)
        , father(-1)
        , leftChild(-1)
        , rightChild(-1)
    {
    }
};

/**
 * @class NewickParser
 * @brief Parses tree structures in Newick format.
 *
 * @details
 * The NewickParser class provides functionality to parse a Newick formatted string
 * and represent the tree structure as a vector of NewickTreeNode structures.
 *
 * @note
 * This parser assumes that the Newick string represents a binary tree. If the tree has nodes
 * with more than two children, the parser will need to be extended accordingly.
 */
class NewickParser {
public:
    NewickParser() = default;
    /**
     * @brief Constructor that takes a Newick formatted string to parse.
     *
     * @param newickStr The Newick format string representing the tree.
     * @throws std::runtime_error If the Newick string has an invalid format.
     */
    explicit NewickParser(const std::string& newickStr);

    void clear() {
        nodes_.clear();
        currentIndex_ = 0;
    }

    /**
     * @brief Retrieves the parsed tree nodes.
     *
     * @return A constant reference to a vector containing all parsed NewickTreeNode structures.
     */
    const std::vector<NewickTreeNode>& getNodes() const;

    std::vector<NewickTreeNode> nodes_;  ///< Vector storing all parsed nodes.
    int currentIndex_ = 0;               ///< Counter to assign unique IDs to nodes.

    /**
     * @brief Main parsing function that initiates the parsing process.
     *
     * @param newickStr The Newick format string to parse.
     * @throws std::runtime_error If the Newick string has an invalid format.
     */
    void parse(const std::string& newickStr);

    /**
     * @brief Recursively parses a subtree in the Newick format.
     *
     * @param str     The Newick string as a character array.
     * @param index   Current parsing position in the string.
     * @param length  Total length of the string.
     * @param father  ID of the parent node. Use -1 if there is no parent (root node).
     * @return The ID of the parsed node.
     * @throws std::runtime_error If the Newick string has an invalid format.
     */
    int parseSubtree(char* str, int& index, int length, int father);

    /**
     * @brief Parses the name of a node from the Newick string.
     *
     * @param str      The Newick string as a character array.
     * @param index    Current parsing position in the string.
     * @param length   Total length of the string.
     * @param outName  Reference to the string where the node name will be stored.
     */
    void parseNodeName(char* str, int& index, int length, std::string& outName);

    /**
     * @brief Parses the branch length (a floating-point number) from the Newick string.
     *
     * @param str      The Newick string as a character array.
     * @param index    Current parsing position in the string.
     * @param length   Total length of the string.
     * @return The parsed branch length as a double.
     * @throws std::runtime_error If the branch length is invalid.
     */
    double parseBranchLength(char* str, int& index, int length);

    /**
     * @brief Skips any whitespace characters in the Newick string.
     *
     * @param str    The Newick string as a character array.
     * @param index  Current parsing position in the string.
     * @param length Total length of the string.
     */
    void skipWhitespace(char* str, int& index, int length);

    /**
     * @brief Trims leading and trailing whitespace from a string.
     *
     * @param s Reference to the string to be trimmed.
     */
    void trimString(std::string& s);

    double distanceBetween(int u, int v) const;

    // NewickParser.hpp 里声明
    std::vector<int> orderLeavesGreedyMinSum(int leafRoot);

    std::vector<std::string> getLeafNames() const;

    int findNodeIdByName(const std::string& name) const;
    void restrictToSubtreeByRootId(int rootId);


};


// -----------------------------
// 工具函数声明
// -----------------------------

// 判断字符串是否是 URL 格式（http/https/ftp）
bool isUrl(const std::string& path_str);

// 验证 URL 是否可以访问
bool verifyUrlReachable(const std::string& url);

// 验证本地文件路径是否存在且合法
void verifyLocalFile(const FilePath& file_path);

// CURL 下载文件时的写入回调函数
size_t writeData(void* ptr, size_t size, size_t nmemb, void* stream);

// 使用 curl 下载远程文件到本地
void downloadFile(const std::string& url, const FilePath& destination);

// 将本地文件复制到目标路径
void copyLocalFile(const FilePath& source, const FilePath& destination);

// 提取文件扩展名字符串
std::string getFileExtension(const FilePath& file_path);

// 拷贝或下载原始数据至工作目录
bool copyRawData(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num);

// 清洗原始数据集，并更新路径映射
bool cleanRawDataset(const FilePath workdir_path, SpeciesPathMap& species_path_map, int thread_num);

// HAL-only preprocessing: keep the existing all-uppercase alignment FASTA and
// record original lowercase runs in a separate export-only sidecar.
bool cleanRawDatasetWithSoftMaskIndex(const FilePath workdir_path,
    SpeciesPathMap& species_path_map,
    SoftMask::PathMap& softmask_path_map,
    int thread_num);

// 运行 WindowMasker 生成重复区域的 interval 文件
std::map<SpeciesName, FilePath> repeatSeqMasking(const FilePath workdir_path, const SpeciesPathMap& species_path_map, int thread_num);

// 根据 interval 文件应用掩码到 FASTA 文件，并更新路径映射
bool applyMaskingAndUpdatePaths(const FilePath workdir_path, SpeciesPathMap& species_path_map, const std::map<SpeciesName, FilePath>& interval_files_map, int thread_num);

// 获取文件大小（支持本地和 URL），输出带单位字符串
std::string getReadableFileSize(const FilePath& filePath);

// 判断文件是否小于指定大小（单位 MB）
bool isFileSmallerThan(const FilePath& filePath, size_t maxSizeMB = 1024);

// 生成一个临时中间处理路径（如 _in_process.fasta）
FilePath getTempFilePath(const FilePath& input_path);

// 获取 .fai 索引的标准路径（原路径加 .fai 后缀）
FilePath getFaiIndexPath(const FilePath& fasta_path);

void repeatMaskRawData(
    const FilePath& work_dir,
    int thread_num,
    SpeciesPathMap& species_path_map
);


// -----------------------------
// seqfile 格式：
// 可选第一行：Newick 进化树字符串
// 其余行：<物种名><空格><FASTA 文件路径 或 URL>
// 例如：
//   ((simGorilla:0.008825,(simHuman:0.0067,simChimp:0.006667)sHuman-sChimp:0.00225)sG-sH-sC:0.00968,simOrang:0.018318);
//   simHuman /path/to/simHuman.fa
//   simChimp http://example.com/simChimp.fa
//   simGorilla /path/to/simGorilla.fa
//   simOrang /path/to/simOrang.fa
//
// 返回 seqfile 中是否包含 Newick 树。无树格式从第一行开始读取物种映射。
// -----------------------------
bool parseSeqfile(const FilePath& seqfile_path,
    NewickParser& newick_tree,
    SpeciesPathMap& species_map,
    const std::string& root = "");

#endif // DATA_PROCESS_H
