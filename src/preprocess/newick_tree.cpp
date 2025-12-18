#include "data_process.h"

#include <vector>
#include <stdexcept>
#include <cctype>    // std::isspace
#include <cstdlib>   // std::strtod
#include <cstring>   // std::strlen
#include <unordered_map>

NewickParser::NewickParser(const std::string& newickStr) {
    // 构造时直接解析 Newick 字符串
    parse(newickStr);
}

const std::vector<NewickTreeNode>& NewickParser::getNodes() const {
    // 返回内部节点数组（只读引用）
    return nodes_;
}

void NewickParser::parse(const std::string& newickStr) {
    // 将字符串复制到可修改的字符缓冲区，便于用指针/下标进行解析
    std::vector<char> buffer(newickStr.begin(), newickStr.end());
    buffer.push_back('\0'); // 末尾补 '\0'，便于 strtod 等函数使用

    int index = 0;
    int length = static_cast<int>(buffer.size()) - 1; // 不包含 '\0' 的有效长度

    // 从根开始解析：father = -1 表示无父节点
    int rootId = parseSubtree(buffer.data(), index, length, -1);

    // 子树解析结束后，跳过空白并尝试匹配分号 ';'
    skipWhitespace(buffer.data(), index, length);
    if (index < length && buffer[index] == ';') {
        index++;
        skipWhitespace(buffer.data(), index, length);
    }

    // 若还有剩余非空白字符，说明 Newick 格式非法
    if (index < length) {
        throw std::runtime_error("Newick format error: Extra characters after parsing.");
    }

    // 若根节点无名称，则默认命名为 "root"
    if (rootId >= 0 && nodes_[rootId].name.empty()) {
        nodes_[rootId].name = "root";
    }
}

int NewickParser::parseSubtree(char* str, int& index, int length, int father) {
    // 跳过空白
    skipWhitespace(str, index, length);

    // 创建一个新节点（先压入占位，后续边解析边填充）
    NewickTreeNode node;
    node.id = currentIndex_++;        // 当前节点编号（单调递增）
    node.father = father;             // 父节点编号
    node.branchLength = 0.0;          // 默认分支长度为 0
    node.isLeaf = false;              // 默认非叶子
    node.leftChild = -1;              // 默认无左孩子
    node.rightChild = -1;             // 默认无右孩子

    // 插入占位节点；后续通过 nodes_[node.id] 回写字段
    nodes_.push_back(node);

    skipWhitespace(str, index, length);

    // ------------------------------
    // 情况 1：内部节点（以 '(' 开头）
    // ------------------------------
    if (index < length && str[index] == '(') {
        index++; // 跳过 '('
        skipWhitespace(str, index, length);

        // 解析左子树
        int leftId = parseSubtree(str, index, length, node.id);
        nodes_[node.id].leftChild = leftId;

        skipWhitespace(str, index, length);

        // 期望遇到 ',' 分隔左右子树
        if (index >= length || str[index] != ',') {
            throw std::runtime_error("Newick format error: Expected ',' after left subtree.");
        }
        index++; // 跳过 ','
        skipWhitespace(str, index, length);

        // 解析右子树
        int rightId = parseSubtree(str, index, length, node.id);
        nodes_[node.id].rightChild = rightId;

        skipWhitespace(str, index, length);

        // 期望遇到 ')' 结束该内部节点的 children 列表
        if (index >= length || str[index] != ')') {
            throw std::runtime_error("Newick format error: Expected ')' after right subtree.");
        }
        index++; // 跳过 ')'
        skipWhitespace(str, index, length);

        // 解析内部节点可选名称
        parseNodeName(str, index, length, nodes_[node.id].name);

        // 解析可选分支长度（以 ':' 开头）
        if (index < length && str[index] == ':') {
            index++; // 跳过 ':'
            skipWhitespace(str, index, length);
            nodes_[node.id].branchLength = parseBranchLength(str, index, length);
        }

        // 明确标记为内部节点
        nodes_[node.id].isLeaf = false;
    }
    // ------------------------------
    // 情况 2：叶节点（不以 '(' 开头）
    // ------------------------------
    else {
        // 标记为叶节点
        nodes_[node.id].isLeaf = true;

        // 解析叶节点名称
        parseNodeName(str, index, length, nodes_[node.id].name);

        // 解析可选分支长度
        if (index < length && str[index] == ':') {
            index++; // 跳过 ':'
            skipWhitespace(str, index, length);
            nodes_[node.id].branchLength = parseBranchLength(str, index, length);
        }

        // 再次确认标记（保持原代码逻辑）
        nodes_[node.id].isLeaf = true;
    }

    // 返回该子树根节点 id
    return node.id;
}

void NewickParser::parseNodeName(char* str, int& index, int length, std::string& outName) {
    // 跳过空白
    skipWhitespace(str, index, length);

    // 节点名称规则：
    // - 允许字母数字与下划线等普通字符
    // - 遇到 ':', ',', ')', '(', ';' 或空白则停止
    while (index < length) {
        char c = str[index];
        if (c == ':' || c == ',' || c == ')' || c == '(' || c == ';' ||
            std::isspace(static_cast<unsigned char>(c))) {
            break;
        }
        outName.push_back(c);
        index++;
    }

    // 额外做一次 trim（保险，避免意外空白）
    trimString(outName);
}

double NewickParser::parseBranchLength(char* str, int& index, int length) {
    // 跳过空白
    skipWhitespace(str, index, length);

    // 记录数字起始指针
    char* startPtr = str + index;
    char* endPtr = nullptr;

    // 使用 strtod 将子串转换为 double
    double val = std::strtod(startPtr, &endPtr);
    if (startPtr == endPtr) {
        // 没有消耗任何字符，说明分支长度非法
        throw std::runtime_error("Newick format error: Invalid branch length.");
    }

    // 推进 index（消耗的字符数 = endPtr - startPtr）
    int consumed = static_cast<int>(endPtr - startPtr);
    index += consumed;

    // 再跳过空白，返回解析到的数值
    skipWhitespace(str, index, length);
    return val;
}

void NewickParser::skipWhitespace(char* str, int& index, int length) {
    // 跳过所有空白字符（空格、制表符、换行等）
    while (index < length && std::isspace(static_cast<unsigned char>(str[index]))) {
        index++;
    }
}

void NewickParser::trimString(std::string& s) {
    if (s.empty()) return;

    // 去掉前导空白
    std::size_t startPos = 0;
    while (startPos < s.size() && std::isspace(static_cast<unsigned char>(s[startPos]))) {
        startPos++;
    }

    // 去掉尾随空白
    std::size_t endPos = s.size();
    while (endPos > startPos && std::isspace(static_cast<unsigned char>(s[endPos - 1]))) {
        endPos--;
    }

    // 截取有效子串
    s = s.substr(startPos, endPos - startPos);
}

// ------------------------------------------------------------------
// 工具函数：返回两个节点之间的距离
// 说明：
// - 使用父指针向上回溯：先记录 u->root 的累计距离，再从 v 向上找 LCA
// - 适合叶子数不大或调用频率不高的场景
// ------------------------------------------------------------------
double NewickParser::distanceBetween(int u, int v) const {
    // 记录从 u 向上到根的路径：distUp[node] = u 到 node 的累积距离
    std::unordered_map<int, double> distUp;
    double acc = 0.0;
    int x = u;
    while (x != -1) {
        distUp[x] = acc;
        int p = nodes_[x].father;
        if (p == -1) break;
        acc += nodes_[x].branchLength; // 注意：沿父边累计长度
        x = p;
    }

    // 再从 v 向上回溯，找到第一个出现在 distUp 的节点，即最近公共祖先（LCA）
    double accV = 0.0;
    int y = v;
    while (y != -1) {
        if (distUp.count(y)) {
            // v->LCA 的距离 + u->LCA 的距离
            return accV + distUp[y];
        }
        int p = nodes_[y].father;
        if (p == -1) break;
        accV += nodes_[y].branchLength;
        y = p;
    }

    // 理论上不会走到这里；作为兜底返回
    return acc + accV;
}

// --------------------------------------------------
// 迭代贪心中心排序：orderLeavesGreedyMinSum
// 目标：
// - 在叶集合上构造一个“从中心到边缘”的顺序
// - 第一个元素是“到所有叶距离之和最小”的叶子（近似树中心）
// 方法：
// - 预计算叶子两两距离矩阵 D
// - 维护每个叶子的 sumDist（到当前集合中其他叶子的距离和）
// - 每次选 sumDist 最小的叶加入 order，并从集合中移除，动态更新其余 sumDist
// --------------------------------------------------
std::vector<int> NewickParser::orderLeavesGreedyMinSum(int rootLeaf) {
    // -------- 1. 收集所有叶节点 --------
    std::vector<int> leaves;
    leaves.reserve(nodes_.size());
    for (const auto& n : nodes_) if (n.isLeaf) leaves.push_back(n.id);

    const int L = static_cast<int>(leaves.size());
    if (L == 0) return {};

    // -------- 2. 预计算 L×L 距离矩阵 --------
    std::vector<std::vector<double>> D(L, std::vector<double>(L, 0.0));
    for (int i = 0; i < L; ++i) {
        for (int j = i + 1; j < L; ++j) {
            double d = distanceBetween(leaves[i], leaves[j]);
            D[i][j] = D[j][i] = d;
        }
    }

    // 叶子 ID -> 距离矩阵索引 的映射
    std::unordered_map<int, int> leaf2idx;
    for (int i = 0; i < L; ++i) leaf2idx[leaves[i]] = i;

    // -------- 3. 初始化 sumDist --------
    // sumDist[i] 表示 leaves[i] 到所有其他叶子的距离之和
    std::vector<double> sumDist(L, 0.0);
    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j)
            if (i != j) sumDist[i] += D[i][j];

    // -------- 4. 迭代贪心选择 --------
    std::vector<int> order;
    order.reserve(L);

    // removed[i] = 1 表示 leaves[i] 已被移除（已选入 order）
    std::vector<char> removed(L, 0);

    for (int step = 0; step < L; ++step) {
        // 在剩余叶中选择 sumDist 最小者
        int bestIdx = -1;
        double bestVal = std::numeric_limits<double>::max();
        for (int i = 0; i < L; ++i) {
            if (removed[i]) continue;
            if (sumDist[i] < bestVal) {
                bestVal = sumDist[i];
                bestIdx = i;
            }
        }

        // 记录该叶子，并将其标记为移除
        order.push_back(leaves[bestIdx]);
        removed[bestIdx] = 1;

        // 更新其余节点的 sumDist：移除 bestIdx 相当于从距离和中减去 D[j][bestIdx]
        for (int j = 0; j < L; ++j) {
            if (removed[j]) continue;
            sumDist[j] -= D[j][bestIdx];
        }
    }

    // 返回顺序：order[0] 是全局距离和最小（中心），后面依次次优
    return order;
}

std::vector<std::string> NewickParser::getLeafNames() const {
    // 返回所有叶节点名称（忽略空名称）
    std::vector<std::string> leaves;
    leaves.reserve(nodes_.size());

    for (const auto& node : nodes_) {
        if (node.isLeaf) {
            if (!node.name.empty()) {
                leaves.push_back(node.name);
            }
        }
    }
    return leaves;
}

// ------------------------------------------------------------------
// 根据节点名称查找节点 id
// - 找到返回 id
// - 未找到返回 -1
// ------------------------------------------------------------------
int NewickParser::findNodeIdByName(const std::string& name) const {
    for (const auto& n : nodes_) {
        if (n.name == name) return n.id;
    }
    return -1;
}

// ------------------------------------------------------------------
// 以指定 rootId 为根裁剪树：仅保留该子树节点，并重建 nodes_ 编号
// 实现要点：
// 1) 从 rootId DFS/BFS 收集所有需保留的 oldId
// 2) oldId -> newId 映射，重建 newNodes（重新编号为 0..k-1）
// 3) 修复 father/leftChild/rightChild 指针到新编号
// 4) 保证新 root 的 father = -1
// ------------------------------------------------------------------
void NewickParser::restrictToSubtreeByRootId(int rootId) {
    if (rootId < 0 || rootId >= (int)nodes_.size())
        throw std::runtime_error("restrictToSubtreeByRootId: invalid root id");

    // -------- 1. 收集子树所有 oldId --------
    std::vector<int> stack{ rootId };
    std::vector<int> kept;
    kept.reserve(nodes_.size());
    std::vector<char> mark(nodes_.size(), 0);

    while (!stack.empty()) {
        int u = stack.back(); stack.pop_back();
        if (mark[u]) continue;
        mark[u] = 1;
        kept.push_back(u);

        const auto& nd = nodes_[u];
        if (nd.leftChild != -1) stack.push_back(nd.leftChild);
        if (nd.rightChild != -1) stack.push_back(nd.rightChild);
    }

    // -------- 2. oldId -> newId 映射，并构建 newNodes --------
    std::vector<int> old2new(nodes_.size(), -1);
    std::vector<NewickTreeNode> newNodes;
    newNodes.reserve(kept.size());

    // 按 kept 顺序重建节点（重新编号）
    for (size_t i = 0; i < kept.size(); ++i) {
        int oldId = kept[i];
        old2new[oldId] = static_cast<int>(i);

        NewickTreeNode nn = nodes_[oldId];
        nn.id = static_cast<int>(i);
        newNodes.push_back(std::move(nn));
    }

    // -------- 3. 修复 parent/children 指针 --------
    for (auto& nn : newNodes) {
        // 修复左/右孩子指针
        nn.leftChild  = (nn.leftChild  != -1 && old2new[nn.leftChild]  != -1) ? old2new[nn.leftChild]  : -1;
        nn.rightChild = (nn.rightChild != -1 && old2new[nn.rightChild] != -1) ? old2new[nn.rightChild] : -1;

        // 修复父指针
        nn.father = (nn.father != -1 && old2new[nn.father] != -1) ? old2new[nn.father] : -1;
    }

    // -------- 4. 确保新 root 的 father == -1 --------
    int newRoot = old2new[rootId];
    if (newRoot >= 0 && newRoot < (int)newNodes.size()) {
        newNodes[newRoot].father = -1;
    }

    // 替换 nodes_ 并更新 currentIndex_
    nodes_.swap(newNodes);
    currentIndex_ = static_cast<int>(nodes_.size());
}
