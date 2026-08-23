#include "align.h"

// ------------------------------------------------------------------
// 函数：cigarToInt
// 功能：将 CIGAR 操作字符（如 'M'）和其长度编码为一个 32 位无符号整数
// 编码格式：高 28 位表示长度，低 4 位表示操作类型编码（参照 SAM 格式）
// 示例：'M', 5 -> 0x00000050 = (5 << 4) | 0x0
// ------------------------------------------------------------------
uint32_t cigarToInt(char operation, uint32_t len) {
	uint32_t opCode;

	// 将字符型操作转换为对应的数字编码（4 位）
	switch (operation) {
	case 'M': opCode = 0x0; break; // 对齐匹配（可能包含错配）
	case 'I': opCode = 0x1; break; // 插入
	case 'D': opCode = 0x2; break; // 删除
	case '=': opCode = 0x7; break; // 精确匹配（无错配）
	case 'X': opCode = 0x8; break; // 错配
		// TODO: 可以根据需要添加更多操作类型（如 soft clip, hard clip 等）
	default: opCode = 0xF; break; // 未知操作，使用保留码 0xF
	}

	// 将长度左移 4 位，然后用 bitwise OR 合并 opCode
	return (len << 4) | opCode;
}

// ------------------------------------------------------------------
// 函数：intToCigar
// 功能：将压缩后的整数值还原为字符型 CIGAR 操作与其长度
// 参数：
//   - cigar：压缩后的整数（低 4 位为 opCode，其他为长度）
//   - operation：输出的操作字符
//   - len：输出的操作长度
// ------------------------------------------------------------------
void intToCigar(CigarUnit cigar, char& operation, uint32_t& len)
{
	uint32_t opCode = cigar & 0xF;  // 提取低 4 位为操作编码
	len = cigar >> 4;               // 高 28 位表示操作长度

	// 将操作编码转换回对应字符
	switch (opCode) {
	case 0x0: operation = 'M'; break;
	case 0x1: operation = 'I'; break;
	case 0x2: operation = 'D'; break;
	case 0x7: operation = '='; break;
	case 0x8: operation = 'X'; break;
	default:  operation = '?'; break; // 未知操作标记为 '?'
	}
}

// -------------------------------------------------------------------
// 函数：caculateMatchScore
// 功能：计算一段对齐序列的总得分（仅基于字符匹配）
// 参数：
//   - match：表示序列的 char 数组（通常来自 CIGAR 'M' 操作对应的比对段）
//   - length：match 的长度
// 返回：整数得分
// 注：此处未使用 HOXD70 替换矩阵，而是使用简单规则：
//     A/T 匹配 +91，C/G 匹配 +100，其他 -100
// ------------------------------------------------------------------
Score_t caculateMatchScore(const char* match, uint_t length)
{
	Score_t score = 0;

	for (size_t i = 0; i < length; ++i) {
		int id = ScoreChar2Idx[match[i]];  // 获取碱基字符在打分矩阵中的索引

		// 根据碱基种类计算得分
		if (id == 0 || id == 3) {            // A 或 T
			score += MATCH_SCORE_AT;
		}
		else if (id == 1 || id == 2) {       // C 或 G
			score += MATCH_SCORE_CG;
		}
		else {                               // 其他如 N
			score += MATCH_SCORE_N;
		}
	}
	return score;
}

/* ------------------------------------------------------------------
 *  追加/拼接 CIGAR：若两端操作码相同则自动合并
 *  ------------------------------------------------------------------
 *  @param dst  目标 CIGAR（被追加）
 *  @param src  待追加的 CIGAR 片段
 * ------------------------------------------------------------------*/
void appendCigar(Cigar_t& dst, const Cigar_t& src)
{
	if (src.empty()) return;                  // nothing to do
	size_t idx = 0;

	//// 1) 检查“拼接点”——dst 的最后一个元素 vs. src 的首元素
	if (!dst.empty()) {
		char op_dst; uint32_t len_dst;
		intToCigar(dst.back(), op_dst, len_dst);

		char op_src; uint32_t len_src;
		intToCigar(src.front(), op_src, len_src);

		// 若操作码相同，则二者合并成一个
		if (op_dst == op_src) {
			dst.back() = cigarToInt(op_dst, len_dst + len_src);
			idx = 1;                         // src 第 0 元素已被合并，后续从 1 开始
		}
	}
	// 2) 将 src 剩余元素逐个 push
	for (; idx < src.size(); ++idx)
		dst.push_back(src[idx]);
}

void prependCigar(Cigar_t& dst, const Cigar_t& src)
{
    if (src.empty()) return; // nothing to do
    size_t end = src.size();

    if (!dst.empty()) {
        char op_dst; uint32_t len_dst;
        intToCigar(dst.front(), op_dst, len_dst);

        char op_src; uint32_t len_src;
        intToCigar(src.back(), op_src, len_src);

        if (op_dst == op_src) {
            // 合并 src.back() 和 dst.front()
            dst.front() = cigarToInt(op_dst, len_dst + len_src);
            end = src.size() - 1; // 最后一个元素已被合并，不再插入
        }
    }

    // 把 src[0 ... end-1] 插到 dst 前面
    dst.insert(dst.begin(), src.begin(), src.begin() + end);
}


/* ------------------------------------------------------------------
 *  追加单个 CIGAR 操作；若与 dst 最后一个操作码一致则合并
 * ------------------------------------------------------------------*/
void appendCigarOp(Cigar_t& dst, char op, uint32_t len)
{
	if (len == 0) return;
	if (!dst.empty()) {
		char op_last; uint32_t len_last;
		intToCigar(dst.back(), op_last, len_last);
		if (op_last == op) {                        // 合并
			dst.back() = cigarToInt(op, len_last + len);
			return;
		}
	}
	dst.push_back(cigarToInt(op, len));
}

/* ----------------------------------------------------------
 *  根据 CIGAR 重建 (ref_aln, qry_aln) 两条同宽序列
 *  --------------------------------------------------------*/
std::pair<std::string, std::string>
buildAlignment(const std::string& ref_seq,
	const std::string& qry_seq,
	const Cigar_t& cigar)
{
	std::string ref_aln;  ref_aln.reserve(ref_seq.size() * 1.2);
	std::string qry_aln;  qry_aln.reserve(qry_seq.size() * 1.2);

	size_t i = 0, j = 0;
	for (CigarUnit u : cigar)
	{
		uint32_t len = u >> 4;
		uint8_t  op = u & 0xF;        // 0=M,1=I,2=D,7='=',8='X'

		if (op == 0 || op == 7 || op == 8) {          // match / mismatch
			ref_aln.append(ref_seq.data() + i, len);
			qry_aln.append(qry_seq.data() + j, len);
			i += len;  j += len;
		}
		else if (op == 1) {                       // insertion in query
			ref_aln.append(len, '-');
			qry_aln.append(qry_seq.data() + j, len);
			j += len;
		}
		else if (op == 2) {                       // deletion in query
			ref_aln.append(ref_seq.data() + i, len);
			qry_aln.append(len, '-');
			i += len;
		}
		// 其余 CIGAR 码如 N/H/S 忽略，本例不出现
	}
	return { ref_aln, qry_aln };
}

/* ----------------------------------------------------------
 *  分割 CIGAR 字符串：根据指定长度将 CIGAR 分为两部分
 *  --------------------------------------------------------*/
std::pair<Cigar_t, Cigar_t> splitCigar(const std::string& cigar_str, uint32_t target_length)
{
    Cigar_t first_part, second_part;
    
    if (cigar_str.empty()) {
        return {first_part, second_part};
    }
    
    if (target_length == 0) {
        // 如果目标长度为0，第一部分为空，第二部分包含整个CIGAR
        size_t pos = 0;
        while (pos < cigar_str.size()) {
            uint32_t op_length = 0;
            while (pos < cigar_str.size() && std::isdigit(cigar_str[pos])) {
                op_length = op_length * 10 + (cigar_str[pos] - '0');
                pos++;
            }
            if (pos >= cigar_str.size()) break;
            char op = cigar_str[pos++];
            second_part.push_back(cigarToInt(op, op_length));
        }
        return {first_part, second_part};
    }
    
    uint32_t current_length = 0;
    size_t pos = 0;
    
    while (pos < cigar_str.size() && current_length < target_length) {
        // 解析数字（长度）
        uint32_t op_length = 0;
        while (pos < cigar_str.size() && std::isdigit(cigar_str[pos])) {
            op_length = op_length * 10 + (cigar_str[pos] - '0');
            pos++;
        }
        
        // 解析操作符
        if (pos >= cigar_str.size()) break;
        char op = cigar_str[pos++];
        
        // 计算此操作在参考序列上消耗的长度
        uint32_t ref_consumed = 0;
        if (op != 'I') {  // 插入操作不消耗参考序列长度
            ref_consumed = op_length;
        }
        
        // 检查是否需要分割这个操作
        if (current_length + ref_consumed <= target_length) {
            // 整个操作都属于第一部分
            first_part.push_back(cigarToInt(op, op_length));
            current_length += ref_consumed;
        } else {
            // 需要分割这个操作
            uint32_t remaining_length = target_length - current_length;
            
            if (op == 'I') {
                // 对于插入操作，如果已经达到目标长度，整个插入都放到第二部分
                second_part.push_back(cigarToInt(op, op_length));
            } else {
                // 对于其他操作，需要分割
                if (remaining_length > 0) {
                    first_part.push_back(cigarToInt(op, remaining_length));
                }
                uint32_t leftover_length = op_length - remaining_length;
                if (leftover_length > 0) {
                    second_part.push_back(cigarToInt(op, leftover_length));
                }
            }
            current_length = target_length;
            break;
        }
    }
    
    // 将剩余的CIGAR操作都加入第二部分
    while (pos < cigar_str.size()) {
        // 解析数字（长度）
        uint32_t op_length = 0;
        while (pos < cigar_str.size() && std::isdigit(cigar_str[pos])) {
            op_length = op_length * 10 + (cigar_str[pos] - '0');
            pos++;
        }
        
        // 解析操作符
        if (pos >= cigar_str.size()) break;
        char op = cigar_str[pos++];
        
        second_part.push_back(cigarToInt(op, op_length));
    }
    
    return {first_part, second_part};
}

/* ----------------------------------------------------------
 *  解析 CIGAR 字符串为 Cigar_t 向量
 *  --------------------------------------------------------*/
void parseCigarString(const std::string& cigar_str, Cigar_t& cigar)
{
    cigar.clear();
    
    if (cigar_str.empty()) return;
    
    size_t pos = 0;
    while (pos < cigar_str.size()) {
        // 解析数字（长度）
        uint32_t op_length = 0;
        while (pos < cigar_str.size() && std::isdigit(cigar_str[pos])) {
            op_length = op_length * 10 + (cigar_str[pos] - '0');
            pos++;
        }
        
        // 解析操作符
        if (pos >= cigar_str.size()) break;
        char op = cigar_str[pos++];
        
        cigar.push_back(cigarToInt(op, op_length));
    }
}

/* ----------------------------------------------------------
 *  将 Cigar_t 转换为字符串格式
 *  --------------------------------------------------------*/
std::string cigarToString(const Cigar_t& cigar)
{
    std::string result;
    for (CigarUnit unit : cigar) {
        char op;
        uint32_t len;
        intToCigar(unit, op, len);
        result += std::to_string(len) + op;
    }
    return result;
}

/* ----------------------------------------------------------
 *  分割 CIGAR 字符串：返回 Cigar_t 和剩余字符串
 *  --------------------------------------------------------*/
std::pair<Cigar_t, std::string> splitCigarMixed(const std::string& cigar_str, uint32_t target_length)
{
    // 先使用原函数获取两个Cigar_t
    auto [first_cigar, second_cigar] = splitCigar(cigar_str, target_length);
    
    // 将第二部分转换为字符串
    std::string second_str = cigarToString(second_cigar);
    
    return {first_cigar, second_str};
}

/* ----------------------------------------------------------
 *  计算 CIGAR 中 M 操作的总长度
 *  --------------------------------------------------------*/
uint32_t countMatchOperations(const Cigar_t& cigar)
{
    uint32_t total_match_length = 0;
    
    for (CigarUnit unit : cigar) {
        char op;
        uint32_t len;
        intToCigar(unit, op, len);
        
        if (op == 'M') {
            total_match_length += len;
        }
    }
    
    return total_match_length;
}

/* ----------------------------------------------------------
 *  计算 CIGAR 中非 D 操作的总长度
 *  --------------------------------------------------------*/
uint32_t countNonDeletionOperations(const Cigar_t& cigar)
{
    uint32_t total_length = 0;
    
    for (CigarUnit unit : cigar) {
        char op;
        uint32_t len;
        intToCigar(unit, op, len);
        
        if (op != 'D') {
            total_length += len;
        }
    }
    
    return total_length;
}

uint32_t countRefLength(const Cigar_t& cigar)
{
    uint32_t total_length = 0;

    for (CigarUnit unit : cigar) {
        char op;
        uint32_t len;
        intToCigar(unit, op, len);

        if (op != 'I') {
            total_length += len;
        }
    }

    return total_length;
}
uint32_t countQryLength(const Cigar_t& cigar)
{
    uint32_t total_length = 0;

    for (CigarUnit unit : cigar) {
        char op;
        uint32_t len;
        intToCigar(unit, op, len);

        if (op != 'D') {
            total_length += len;
        }
    }

    return total_length;
}


uint32_t countAlignmentLength(const Cigar_t& cigar) {
    uint32_t total = 0;
    for (auto unit : cigar) {
        uint32_t len = unit >> 4;   // 高 28 位是长度
        total += len;
    }
    return total;
}


/// \brief 检查 gap_cigar 的质量
/// \param cigar   WFA/KSW2 生成的 Cigar_t
/// \param ref_len gap 区间 ref 的长度
/// \param qry_len gap 区间 query 的长度
/// \return 是否达标（true=通过，false=不通过）
bool checkGapCigarQuality(const Cigar_t& cigar,
    size_t ref_len,
    size_t qry_len,
    double min_identity) {
    if (cigar.empty()) return false;

    size_t max_gap_run = 0;    // 当前连续 gap 的长度
    size_t matches = 0;    // 匹配长度
    size_t aln_len = 0;    // 对齐长度（不计 clipping）

    for (auto unit : cigar) {
        uint32_t len;
        char op;
        intToCigar(unit, op, len);

        if (op == 'M' || op == '=' || op == 'X') {
            // 匹配/错配
            aln_len += len;
            if (op == 'M' || op == '=') matches += len;
            max_gap_run = 0; // 断开 gap run
        }
        else if (op == 'I' || op == 'D') {
            aln_len += len;
            max_gap_run += len;

        }
        else {
            // 其他操作（S, H, N等），这里直接跳过
            continue;
        }
    }

    if (aln_len <= 10) return true;
	uint_t min_len = std::min(ref_len, qry_len);
    double identity = static_cast<double>(matches) / min_len;

    return (identity >= min_identity);
}

bool alignmentRowsPreferredToUnaligned(const std::string& reference,
                                       const std::string& query) {
    if (reference.empty() || reference.size() != query.size()) {
        return false;
    }

    enum class GapSide : uint8_t {
        NONE = 0,
        REFERENCE,
        QUERY
    };
    GapSide gap_side = GapSide::NONE;
    int64_t alignment_score = 0;
    size_t reference_bases = 0;
    size_t query_bases = 0;
    for (size_t index = 0; index < reference.size(); ++index) {
        const bool reference_gap = reference[index] == '-';
        const bool query_gap = query[index] == '-';
        if (reference_gap && query_gap) {
            continue;
        }
        if (reference_gap || query_gap) {
            const GapSide current_side =
                reference_gap ? GapSide::REFERENCE : GapSide::QUERY;
            alignment_score -= current_side == gap_side
                                   ? GAP_EXTEND_PENALTY
                                   : GAP_OPEN_PENALTY;
            gap_side = current_side;
            reference_bases += !reference_gap;
            query_bases += !query_gap;
            continue;
        }
        gap_side = GapSide::NONE;
        alignment_score += subsScore(reference[index], query[index]);
        ++reference_bases;
        ++query_bases;
    }

    const auto unaligned_gap_score = [](size_t length) -> int64_t {
        if (length == 0) {
            return 0;
        }
        return -static_cast<int64_t>(GAP_OPEN_PENALTY) -
               static_cast<int64_t>(length - 1) *
                   GAP_EXTEND_PENALTY;
    };
    const int64_t unaligned_score =
        unaligned_gap_score(reference_bases) +
        unaligned_gap_score(query_bases);
    return alignment_score > unaligned_score;
}

bool alignmentCigarPreferredToUnaligned(
    const std::string& reference,
    const std::string& query,
    const Cigar_t& cigar) {
    if (cigar.empty()) {
        return false;
    }

    enum class GapSide : uint8_t {
        NONE = 0,
        REFERENCE,
        QUERY
    };
    const auto unaligned_gap_score = [](size_t length) -> int64_t {
        if (length == 0) {
            return 0;
        }
        return -static_cast<int64_t>(GAP_OPEN_PENALTY) -
               static_cast<int64_t>(length - 1) *
                   GAP_EXTEND_PENALTY;
    };

    GapSide gap_side = GapSide::NONE;
    int64_t alignment_score = 0;
    size_t reference_offset = 0;
    size_t query_offset = 0;
    size_t alignment_columns = 0;
    for (const CigarUnit unit : cigar) {
        const uint32_t length = unit >> 4;
        const uint8_t operation = unit & 0xF;
        if (length == 0) {
            continue;
        }
        if (operation == 0 || operation == 7 || operation == 8) {
            if (length > reference.size() - reference_offset ||
                length > query.size() - query_offset) {
                return false;
            }
            gap_side = GapSide::NONE;
            for (uint32_t index = 0; index < length; ++index) {
                alignment_score += subsScore(
                    reference[reference_offset + index],
                    query[query_offset + index]);
            }
            reference_offset += length;
            query_offset += length;
            alignment_columns += length;
        } else if (operation == 1) {
            if (length > query.size() - query_offset) {
                return false;
            }
            alignment_score -=
                gap_side == GapSide::REFERENCE
                    ? static_cast<int64_t>(length) * GAP_EXTEND_PENALTY
                    : static_cast<int64_t>(GAP_OPEN_PENALTY) +
                          static_cast<int64_t>(length - 1) *
                              GAP_EXTEND_PENALTY;
            gap_side = GapSide::REFERENCE;
            query_offset += length;
            alignment_columns += length;
        } else if (operation == 2) {
            if (length > reference.size() - reference_offset) {
                return false;
            }
            alignment_score -=
                gap_side == GapSide::QUERY
                    ? static_cast<int64_t>(length) * GAP_EXTEND_PENALTY
                    : static_cast<int64_t>(GAP_OPEN_PENALTY) +
                          static_cast<int64_t>(length - 1) *
                              GAP_EXTEND_PENALTY;
            gap_side = GapSide::QUERY;
            reference_offset += length;
            alignment_columns += length;
        }
    }

    if (alignment_columns == 0) {
        return false;
    }
    const int64_t unaligned_score =
        unaligned_gap_score(reference_offset) +
        unaligned_gap_score(query_offset);
    return alignment_score > unaligned_score;
}

