// =============================================================
//  File: ramesh.cpp   –  Low‑level primitives (v0.6‑alpha)
//  ✧  Segment / Block / GenomeEnd implementation  ✧
// =============================================================
#include "extension_internal.h"
#include "ramesh.h"
#include <algorithm>
#include <iomanip>
#include <shared_mutex>
#include <stdexcept>
#include <type_traits>

namespace RaMesh {
    namespace {
        std::atomic<uint64_t> g_next_block_id{1};

        Segment* findExtensionCandidate(
            Segment* reference_segment,
            const SpeciesName& query_name,
            bool is_left_extend) {
            if (!reference_segment) return nullptr;

            SegPtr candidate = is_left_extend
                ? reference_segment->primary_path.prev.load(
                    std::memory_order_acquire)
                : reference_segment->primary_path.next.load(
                    std::memory_order_acquire);
            while (candidate) {
                if ((is_left_extend && candidate->isHead()) ||
                    (!is_left_extend && candidate->isTail())) {
                    break;
                }
                // Preserve the legacy ordering: a fully extended segment is a
                // barrier even if its block also contains query_name.
                if (candidate->left_extend && candidate->right_extend) {
                    break;
                }

                const Block* block = candidate->parent_block.get();
                if (!block) {
                    throw std::runtime_error(
                        "Reference candidate is missing its parent block");
                }
                for (const auto& [key, segment] : block->anchors) {
                    (void)segment;
                    if (key.first == query_name) {
                        return candidate.get();
                    }
                }
                candidate = is_left_extend
                    ? candidate->primary_path.prev.load(
                        std::memory_order_acquire)
                    : candidate->primary_path.next.load(
                        std::memory_order_acquire);
            }
            return nullptr;
        }
    }

    namespace detail {

    const SeqPro::SequenceManager& originalSequenceManager(
        const SeqPro::SharedManagerVariant& manager) {
        if (!manager) {
            throw std::invalid_argument("Sequence manager is null");
        }
        return std::visit(
            [](const auto& pointer) -> const SeqPro::SequenceManager& {
                if (!pointer) {
                    throw std::invalid_argument(
                        "Sequence manager variant contains a null pointer");
                }
                using Pointer = std::decay_t<decltype(pointer)>;
                if constexpr (std::is_same_v<
                                  Pointer,
                                  std::unique_ptr<SeqPro::SequenceManager>>) {
                    return *pointer;
                } else {
                    return pointer->getOriginalManager();
                }
            },
            *manager);
    }

    PreparedExtensionResult alignIntervalPrepared(
        Segment* query_segment,
        Segment* reference_segment,
        Segment* reference_candidate,
        const ExtensionSequenceContext& query_context,
        const ExtensionSequenceContext& reference_context,
        bool is_left_extend,
        int_t zdrop,
        ExtensionSequenceScratch& scratch) {
        if (!query_segment || !reference_segment ||
            query_segment->isHead() || query_segment->isTail() ||
            (query_segment->left_extend && query_segment->right_extend) ||
            (is_left_extend && query_segment->left_extend) ||
            (!is_left_extend && query_segment->right_extend)) {
            return {PreparedExtensionStatus::SKIPPED};
        }
        if (!reference_candidate) {
            return {PreparedExtensionStatus::NO_CANDIDATE};
        }

        int_t query_start = 0;
        int_t query_length = 0;
        int_t reference_start = 0;
        int_t reference_length = 0;
        const Strand strand = query_segment->strand;

        if (is_left_extend) {
            if (strand == FORWARD) {
                const SegPtr query_left =
                    query_segment->primary_path.prev.load(
                        std::memory_order_acquire);
                query_start = !query_left->isHead()
                    ? static_cast<int_t>(
                          query_left->start + query_left->length)
                    : 0;
                query_length = static_cast<int_t>(query_segment->start) -
                    query_start;
            } else {
                const SegPtr query_right =
                    query_segment->primary_path.next.load(
                        std::memory_order_acquire);
                query_start = static_cast<int_t>(
                    query_segment->start + query_segment->length);
                query_length = !query_right->isTail()
                    ? static_cast<int_t>(query_right->start) - query_start
                    : static_cast<int_t>(query_context.chromosome_length) -
                          query_start;
            }

            reference_start = static_cast<int_t>(
                reference_candidate->start + reference_candidate->length);
            reference_length = static_cast<int_t>(reference_segment->start) -
                reference_start;
        } else {
            if (strand == FORWARD) {
                const SegPtr query_right =
                    query_segment->primary_path.next.load(
                        std::memory_order_acquire);
                query_start = static_cast<int_t>(
                    query_segment->start + query_segment->length);
                query_length = !query_right->isTail()
                    ? static_cast<int_t>(query_right->start) - query_start
                    : static_cast<int_t>(query_context.chromosome_length) -
                          query_start;
            } else {
                const SegPtr query_left =
                    query_segment->primary_path.prev.load(
                        std::memory_order_acquire);
                query_start = !query_left->isHead()
                    ? static_cast<int_t>(
                          query_left->start + query_left->length)
                    : 0;
                query_length = static_cast<int_t>(query_segment->start) -
                    query_start;
            }

            reference_start = static_cast<int_t>(
                reference_segment->start + reference_segment->length);
            reference_length =
                static_cast<int_t>(reference_candidate->start) -
                reference_start;
        }

        if (query_length <= 0 || reference_length <= 0) {
            return {PreparedExtensionStatus::NONPOSITIVE_GAP};
        }
        if (query_length > kMaximumExtensionGap ||
            reference_length > kMaximumExtensionGap) {
            return {PreparedExtensionStatus::GAP_TOO_LONG};
        }
        if (!query_context.manager || !reference_context.manager) {
            throw std::invalid_argument(
                "Prepared extension is missing a sequence manager");
        }

        query_context.manager->getSubSequenceInto(
            query_context.sequence_id,
            static_cast<SeqPro::Position>(
                static_cast<Coord_t>(query_start)),
            static_cast<SeqPro::Length>(
                static_cast<Coord_t>(query_length)),
            scratch.query);
        reference_context.manager->getSubSequenceInto(
            reference_context.sequence_id,
            static_cast<SeqPro::Position>(
                static_cast<Coord_t>(reference_start)),
            static_cast<SeqPro::Length>(
                static_cast<Coord_t>(reference_length)),
            scratch.reference);

        if (is_left_extend) {
            std::reverse(scratch.reference.begin(), scratch.reference.end());
            if (strand == FORWARD) {
                std::reverse(scratch.query.begin(), scratch.query.end());
            } else {
                baseComplement(scratch.query);
            }
        } else if (strand == REVERSE) {
            reverseComplement(scratch.query);
        }

        Cigar_t result =
            extendAlignKSW2(scratch.reference, scratch.query, zdrop);
        if (!alignmentCigarPreferredToUnaligned(
                scratch.reference, scratch.query, result)) {
            return {PreparedExtensionStatus::ALIGNMENT_REJECTED};
        }

        if (is_left_extend) {
            std::reverse(result.begin(), result.end());
            const AlignCount count = countAlignedBases(result);
            if (strand == FORWARD) {
                query_segment->start -= count.query_bases;
            }
            query_segment->length += count.query_bases;
            reference_segment->start -= count.ref_bases;
            reference_segment->length += count.ref_bases;
            prependCigar(query_segment->cigar, result);
        } else {
            const AlignCount count = countAlignedBases(result);
            if (strand == REVERSE) {
                query_segment->start -= count.query_bases;
            }
            query_segment->length += count.query_bases;
            reference_segment->length += count.ref_bases;
            appendCigar(query_segment->cigar, result);
        }
        return {PreparedExtensionStatus::ALIGNMENT_ACCEPTED};
    }

    }  // namespace detail


    /* =============================================================
     * 0.  Segment factories & utilities
     * ===========================================================*/
    SegPtr Segment::create(uint_t start, uint_t len, Strand sd,
        Cigar_t cg, AlignRole rl, SegmentRole sl,
        const BlockPtr& bp)
    {
        auto s = std::make_shared<Segment>();
        s->start = start;
        s->length = len;
        s->strand = sd;
        s->cigar = std::move(cg);
        s->align_role = rl;
        s->seg_role = sl;
        if (bp) s->parent_block = bp;

        // list pointers – initialise nullptr
        s->primary_path.next.store(nullptr, std::memory_order_relaxed);
        s->primary_path.prev.store(nullptr, std::memory_order_relaxed);

        return s;
    }

    SegPtr Segment::createFromRegion(Region& region, Strand sd,
        Cigar_t cg, AlignRole rl, SegmentRole sl,
        const BlockPtr& bp)
    {
        return create(region.start, region.length, sd, std::move(cg), rl, sl, bp);
    }

    SegPtr Segment::createHead()
    {
        auto* h = new Segment();
        h->seg_role = SegmentRole::HEAD;
        h->align_role = AlignRole::PRIMARY;
        h->primary_path.next.store(nullptr, std::memory_order_relaxed);
        h->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return std::shared_ptr<Segment>(h);
    }

    SegPtr Segment::createTail()
    {
        auto* t = new Segment();
        t->seg_role = SegmentRole::TAIL;
        t->align_role = AlignRole::PRIMARY;
        t->primary_path.next.store(nullptr, std::memory_order_relaxed);
        t->primary_path.prev.store(nullptr, std::memory_order_relaxed);
        return std::shared_ptr<Segment>(t);
    }

    //void Segment::linkChain(const std::vector<SegPtr>& segs)
    //{
    //    for (size_t i = 0; i + 1 < segs.size(); ++i) {
    //        segs[i]->primary_path.next.store(segs[i + 1], std::memory_order_relaxed);
    //        segs[i + 1]->primary_path.prev.store(segs[i], std::memory_order_release);
    //    }
    //}
    void Segment::linkChain(const std::vector<SegPtr>& segs)
    {
        for (size_t i = 0; i + 1 < segs.size(); ++i) {
            segs[i]->primary_path.next.store(segs[i + 1], std::memory_order_release);
            segs[i + 1]->primary_path.prev.store(segs[i], std::memory_order_release);
        }
    }

    
    //void Segment::unlinkSegment(SegPtr segment) {
    //    if (!segment || segment->isHead() || segment->isTail()) return;
    //    
    //    // 获取前驱和后继
    //    SegPtr prev = segment->primary_path.prev.load(std::memory_order_acquire);
    //    SegPtr next = segment->primary_path.next.load(std::memory_order_acquire);
    //    
    //    if (prev && next) {
    //        // 原子地更新链表指针
    //        prev->primary_path.next.store(next, std::memory_order_release);
    //        next->primary_path.prev.store(prev, std::memory_order_release);
    //        
    //        // 清空被删除segment的指针
    //        segment->primary_path.next.store(nullptr, std::memory_order_relaxed);
    //        segment->primary_path.prev.store(nullptr, std::memory_order_relaxed);
    //    }
    //}
// 把  | prev ⇆ seg ⇆ next |  改成  | prev ⇆ next |
    void Segment::unlinkSegment(SegPtr seg)
    {
        if (!seg || seg->isHead() || seg->isTail()) return;

        SegPtr prev = seg->primary_path.prev.load(std::memory_order_acquire);
        SegPtr next = seg->primary_path.next.load(std::memory_order_acquire);

        if (!prev || !next) {
            seg->primary_path.next.store(nullptr, std::memory_order_release);
            seg->primary_path.prev.store(nullptr, std::memory_order_release);
            return;
        }

        if (prev->primary_path.next.load(std::memory_order_acquire) == seg &&
            next->primary_path.prev.load(std::memory_order_acquire) == seg) {
            prev->primary_path.next.store(next, std::memory_order_release);
            next->primary_path.prev.store(prev, std::memory_order_release);
        }

        seg->primary_path.next.store(nullptr, std::memory_order_release);
        seg->primary_path.prev.store(nullptr, std::memory_order_release);
    }


    
    void Segment::deleteSegment(SegPtr segment) {
        if (!segment) return;
        
        // 先从链表中解除链接
        unlinkSegment(segment);
        
        // 清理parent_block引用
        if (segment->parent_block) {
            std::unique_lock block_lock(segment->parent_block->rw);
            
            // 从block的anchors中移除此segment
            for (auto it = segment->parent_block->anchors.begin(); 
                 it != segment->parent_block->anchors.end(); ++it) {
                if (it->second == segment) {
                    segment->parent_block->anchors.erase(it);
                    break;
                }
            }
            segment->parent_block.reset();
        }
        
    }
    
    void Segment::deleteBatch(const std::vector<SegPtr>& segments) {
        if (segments.empty()) return;
        
        // 按block分组以减少锁竞争
        std::unordered_map<BlockPtr, std::vector<SegPtr>> block_groups;
        std::vector<SegPtr> orphaned_segments;
        
        for (SegPtr seg : segments) {
            if (!seg || seg->isHead() || seg->isTail()) continue;
            
            if (seg->parent_block) {
                block_groups[seg->parent_block].emplace_back(seg);
            } else {
                orphaned_segments.emplace_back(seg);
            }
        }
        
        // 批量处理每个block
        for (auto& [block, segs] : block_groups) {
            std::unique_lock block_lock(block->rw);
            
            for (SegPtr seg : segs) {
                // 从链表中解除链接
                unlinkSegment(seg);
                
                // 从anchors中移除
                for (auto it = block->anchors.begin(); it != block->anchors.end(); ++it) {
                    if (it->second == seg) {
                        block->anchors.erase(it);
                        break;
                    }
                }
                seg->parent_block.reset();
            }
        }
        
        // 处理孤立的segments
        for (SegPtr seg : orphaned_segments) {
            unlinkSegment(seg);
        }
    }

    /* =============================================================
     * 1. GenomeEnd helpers
     * ===========================================================*/
    GenomeEnd::GenomeEnd() {
        head = Segment::createHead();
        tail = Segment::createTail();
        head->primary_path.next.store(tail, std::memory_order_relaxed);
        tail->primary_path.prev.store(head, std::memory_order_release);

        sample_vec.resize(1, head);   // slot 0 永远指向 head
    }

    /* ---------- 采样表维护 ---------- */
    void GenomeEnd::ensureSampleSize(uint_t pos) {
        std::size_t need = pos / kSampleStep + 1;
        if (need > sample_vec.size()) sample_vec.resize(need, nullptr);
    }

    void GenomeEnd::setToSampling(SegPtr cur) {
        // std::unique_lock lk(rw);

        std::size_t idx = cur->start / kSampleStep;
        
        if (idx == 0) {
            return;
        }
        if (idx + 1 > sample_vec.size()) sample_vec.resize(idx + 1, nullptr);

        if (!sample_vec[idx] || !sample_vec[idx]->parent_block || cur->start > sample_vec[idx]->start)
        {
            sample_vec[idx] = cur;
        }

    }

    void GenomeEnd::updateSampling(const std::vector<SegPtr>& segs) {
        if (segs.empty()) return;
        std::unique_lock lk(rw);      // ××× 写锁 —— 修改 sample_vec
        for (SegPtr s : segs) {
            std::size_t idx = s->start / kSampleStep;
            ensureSampleSize(s->start);
            if (!sample_vec[idx] || s->start < sample_vec[idx]->start)
                sample_vec[idx] = s;   // 只保留区间内最左端 segment
        }
    }

    SegPtr GenomeEnd::findSurrounding(uint_t range_start) {
        // 1) 读取采样表得到“最近前驱”的 hint
        //std::shared_lock lk(rw);                 // 读锁即可
        std::size_t slot = std::max((size_t)(range_start / kSampleStep) - 1, (size_t)0);
        //lk.unlock();                             // 之后只读链表，不再访问 sample_vec
        SegPtr hint = (slot < sample_vec.size() && sample_vec[slot])
            ? sample_vec[slot] : head;
   
        // 2) 保证 hint 在目标区间左侧
        while (!hint->isHead() && hint->start > range_start) {
            SegPtr nxt = hint->primary_path.prev.load(std::memory_order_acquire);
            if (!nxt) {                    // 保险：有人在并发删除
                hint = head;               // 回到链表起点重新来
                break;
            }
            hint = nxt;
        }

        SegPtr prev = hint;
        SegPtr curr = hint->primary_path.next.load(std::memory_order_acquire);

        // 3) 向右遍历，直到越过 range_start
        while (curr && !curr->isTail() && curr->start <= range_start) {
            prev = curr;
            curr = curr->primary_path.next.load(std::memory_order_acquire);
        }
        return prev;
    }




    void GenomeEnd::insertSegment(const SegPtr seg)
    {
        if (!seg) return;

        uint_t beg = seg->start;

        // 1) 找到目标区间的前驱/后继（只读操作）
        SegPtr prev = findSurrounding(beg);
        SegPtr next = prev->primary_path.next.load();

        seg->primary_path.prev.store(prev);
        
        seg->primary_path.next.store(next);
      
        prev->primary_path.next.store(seg);
        next->primary_path.prev.store(seg);
       
        setToSampling(seg);
    }

    void GenomeEnd::insertSegmentsSorted(
        const std::vector<SegPtr>& segments)
    {
        if (segments.empty()) return;
        SegPtr previous;
        uint_t previous_start = 0;
        bool first = true;
        for (const SegPtr& segment : segments) {
            if (!segment) continue;
            if (first || segment->start - previous_start > 2 * kSampleStep) {
                previous = findSurrounding(segment->start);
            } else {
                SegPtr current = previous->primary_path.next.load(
                    std::memory_order_acquire);
                while (current && !current->isTail() &&
                       current->start <= segment->start) {
                    previous = current;
                    current = current->primary_path.next.load(
                        std::memory_order_acquire);
                }
            }
            SegPtr next = previous->primary_path.next.load(
                std::memory_order_acquire);
            segment->primary_path.prev.store(
                previous, std::memory_order_release);
            segment->primary_path.next.store(
                next, std::memory_order_release);
            previous->primary_path.next.store(
                segment, std::memory_order_release);
            next->primary_path.prev.store(
                segment, std::memory_order_release);
            setToSampling(segment);
            previous = segment;
            previous_start = segment->start;
            first = false;
        }
    }


    void GenomeEnd::clearAllSegments() {
        std::unique_lock lk(rw);

        // 1) 复位双向链表
        head->primary_path.next.store(tail, std::memory_order_relaxed);
        tail->primary_path.prev.store(head, std::memory_order_release);

        // 2) 清空采样表，只保留 slot0 指向 head
        sample_vec.clear();
        sample_vec.resize(1, head);
    }
    
    bool GenomeEnd::removeSegment(SegPtr segment) {
        if (!segment || segment->isHead() || segment->isTail()) return false;
        
        std::unique_lock lk(rw);
        
        // 验证segment确实在这个链表中
        SegPtr current = head->primary_path.next.load(std::memory_order_acquire);
        bool found = false;
        while (current && !current->isTail()) {
            if (current == segment) {
                found = true;
                break;
            }
            current = current->primary_path.next.load(std::memory_order_acquire);
        }
        
        if (!found) return false;
        
        // 从链表中移除
        Segment::unlinkSegment(segment);
        
        // 更新采样表
        invalidateSampling(segment->start, segment->start + segment->length);
        
        return true;
    }
    
    bool GenomeEnd::removeSegmentRange(uint_t range_start, uint_t range_end) {
        std::unique_lock lk(rw);
        
        std::vector<SegPtr> to_remove;
        SegPtr current = head->primary_path.next.load(std::memory_order_acquire);
        
        // 收集需要删除的segments
        while (current && !current->isTail()) {
            if (current->start >= range_start && 
                current->start + current->length <= range_end) {
                to_remove.emplace_back(current);
            }
            current = current->primary_path.next.load(std::memory_order_acquire);
        }
        
        // 批量删除
        for (SegPtr seg : to_remove) {
            Segment::unlinkSegment(seg);
        }
        
        if (!to_remove.empty()) {
            invalidateSampling(range_start, range_end);
        }
        
        return !to_remove.empty();
    }
    
    void GenomeEnd::removeBatch(const std::vector<SegPtr>& segments) {
        if (segments.empty()) return;
        
        std::unique_lock lk(rw);
        
        uint_t min_pos = UINT32_MAX;
        uint_t max_pos = 0;
        
        // 批量删除segments并记录范围
        for (SegPtr seg : segments) {
            if (seg && !seg->isHead() && !seg->isTail()) {
                Segment::unlinkSegment(seg);
                min_pos = std::min(min_pos, seg->start);
                max_pos = std::max(max_pos, seg->start + seg->length);
            }
        }
        
        // 更新采样表
        if (min_pos != UINT32_MAX) {
            invalidateSampling(min_pos, max_pos);
        }
    }
    
    void GenomeEnd::invalidateSampling(uint_t start, uint_t end) {
        // 计算受影响的采样区间
        std::size_t start_idx = start / kSampleStep;
        std::size_t end_idx = end / kSampleStep + 1;
        
        // 重建受影响区间的采样
        for (std::size_t idx = start_idx; idx <= end_idx && idx < sample_vec.size(); ++idx) {
            sample_vec[idx] = nullptr;
            
            // 重新扫描该区间找到最左的segment
            uint_t region_start = idx * kSampleStep;
            uint_t region_end = (idx + 1) * kSampleStep;
            
            SegPtr current = head->primary_path.next.load(std::memory_order_acquire);
            while (current && !current->isTail()) {
                if (current->start >= region_start && current->start < region_end) {
                    if (!sample_vec[idx] || current->start < sample_vec[idx]->start) {
                        sample_vec[idx] = current;
                    }
                }
                current = current->primary_path.next.load(std::memory_order_acquire);
            }
        }
    }

    void RaMesh::GenomeEnd::alignInterval(
        const SpeciesName& ref_name,
        const SpeciesName& query_name,
        const ChrName& query_chr_name,
        SegPtr cur_node,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& managers,
        bool is_left_extend,
        int_t zdrop) {
        if (!cur_node || cur_node == head || cur_node == tail ||
            (cur_node->right_extend && cur_node->left_extend) ||
            (is_left_extend && cur_node->left_extend) ||
            (!is_left_extend && cur_node->right_extend)) {
            return;
        }

        Block* current_block = cur_node->parent_block.get();
        if (!current_block) {
            throw std::runtime_error(
                "Query segment is missing its parent block");
        }
        const auto reference_occurrence = current_block->anchors.find(
            {ref_name, current_block->ref_chr});
        if (reference_occurrence == current_block->anchors.end()) {
            throw std::runtime_error(
                "Reference occurrence is missing from the current block");
        }

        Segment* reference_segment = reference_occurrence->second.get();
        Segment* reference_candidate = findExtensionCandidate(
            reference_segment, query_name, is_left_extend);
        if (!reference_candidate) return;

        const auto& query_manager = detail::originalSequenceManager(
            managers.at(query_name));
        const auto& reference_manager = detail::originalSequenceManager(
            managers.at(ref_name));
        const SeqPro::SequenceId query_sequence_id =
            query_manager.getSequenceId(query_chr_name);
        const SeqPro::SequenceId reference_sequence_id =
            reference_manager.getSequenceId(current_block->ref_chr);

        const detail::ExtensionSequenceContext query_context{
            &query_manager,
            query_sequence_id,
            query_manager.getSequenceLength(query_sequence_id)};
        const detail::ExtensionSequenceContext reference_context{
            &reference_manager,
            reference_sequence_id,
            reference_manager.getSequenceLength(reference_sequence_id)};
        detail::ExtensionSequenceScratch scratch;
        scratch.reserve();
        (void)detail::alignIntervalPrepared(
            cur_node.get(), reference_segment, reference_candidate,
            query_context, reference_context, is_left_extend, zdrop, scratch);
    }

    // 串行重新排序：按照 start 坐标
    void GenomeEnd::resortSegments() {
        SegPtr cur = head->primary_path.next.load(std::memory_order_acquire);
        if (!cur || cur->isTail()) return;

        uint_t previous_start = cur->start;
        cur = cur->primary_path.next.load(std::memory_order_acquire);
        bool strictly_increasing = true;
        while (cur && !cur->isTail()) {
            if (cur->start <= previous_start) {
                strictly_increasing = false;
                break;
            }
            previous_start = cur->start;
            cur = cur->primary_path.next.load(std::memory_order_acquire);
        }
        if (strictly_increasing && cur && cur->isTail()) {
            return;
        }

        std::vector<SegPtr> segs;
        cur = head->primary_path.next.load(std::memory_order_acquire);
        while (cur && !cur->isTail()) {
            segs.push_back(cur);
            cur = cur->primary_path.next.load(std::memory_order_acquire);
        }

        if (segs.empty()) return;

        // 按 start 排序
        std::sort(segs.begin(), segs.end(), [](const SegPtr& a, const SegPtr& b) {
            return a->start < b->start;
            });

        // 重新链接：head -> segs[0] -> ... -> segs[n-1] -> tail
        Segment::linkChain(segs);
        head->primary_path.next.store(segs.front(), std::memory_order_release);
        segs.front()->primary_path.prev.store(head, std::memory_order_release);

        tail->primary_path.prev.store(segs.back(), std::memory_order_release);
        segs.back()->primary_path.next.store(tail, std::memory_order_release);

#ifdef _DEBUG_
        cur = head->primary_path.next.load(std::memory_order_acquire);
        uint_t prev_start = 0;
        bool first = true;
        while (cur && !cur->isTail()) {
            if (!first && cur->start < prev_start) {
                std::cerr << "[GenomeEnd::resortSegments] ERROR: "
                    << "segment order invalid: "
                    << cur->start << " < " << prev_start << std::endl;
                break; // 发现问题就退出
            }
            prev_start = cur->start;
            first = false;
            cur = cur->primary_path.next.load(std::memory_order_acquire);
        }
#endif
    }


    /* -------------------------------------------------------------
     *  移除同一染色体链表中的重叠 Segment
     *  规则：两条 Segment 区间有任何交叠时，删除 length 较小者
     * -------------------------------------------------------------*/
    void RaMesh::GenomeEnd::removeOverlap(bool if_ref)
    {
        std::unique_lock<std::shared_mutex> lk(rw);      // 独占写


        if (!head || !tail) return;

        SegPtr prev = head->primary_path.next.load(std::memory_order_acquire);
        if (!prev || prev == tail) return;

        SegPtr cur = prev->primary_path.next.load(std::memory_order_acquire);

        while (cur && cur != tail)
        {
            SegPtr next = cur;
            uint_t prev_end = prev->start + prev->length;
            uint_t cur_start = cur->start;
 
            /* -------- 检测交叠 -------- */
            while (prev_end > cur_start && cur && cur != tail && next && next != tail)
            {
                cur_start = cur->start;

                if (if_ref && prev->cigar.size() == 0 && cur->cigar.size() == 0) {
                    cur = cur->primary_path.next.load(std::memory_order_acquire);
                }
                else {
                    bool cur_longer = prev->length <= cur->length;
                    if (cur_longer) {
                        prev->parent_block->removeAllSegments();
                        break;
                    }
                    else {
                        if (next == cur) {
                            next = cur->primary_path.next.load(std::memory_order_acquire);
                        }
                        cur = cur->primary_path.next.load(std::memory_order_acquire);                   
                        SegPtr cur_prev = cur->primary_path.prev.load(std::memory_order_acquire);                       
                        cur_prev->parent_block->removeAllSegments();
                    }
                }
            }
            // setToSampling(next);
            /* -------- 无交叠，正常前进 -------- */
            prev = next;
            cur = prev->primary_path.next.load(std::memory_order_acquire);
        }

    }
    /* -------------------------------------------------------------
 *  移除同一染色体链表中的重叠 Segment
 *  规则：两条 Segment 区间有任何交叠时，删除 length 较小者
 * ------------------------------------------------------------*/
    //void GenomeEnd::removeOverlap(bool if_ref)
    //{
    //    std::unique_lock lk(rw);                 // 串行调用，独占即可

    //    SegPtr seg = head->primary_path.next.load(std::memory_order_relaxed);
    //    while (seg && !seg->isTail())
    //    {
    //        SegPtr nxt = seg->primary_path.next.load(std::memory_order_relaxed);

    //        /* 检查 seg 与后继是否重叠 —— 只要起点落在 seg 区间内就算重叠 */
    //        while (nxt && !nxt->isTail() &&
    //            nxt->start < seg->start + seg->length)
    //        {
    //            /* 如果是 ref 链且两段都是“纯 M”（cigar 为空），
    //               你曾选择保留两段，这里沿用原逻辑                 */
    //            if (if_ref && seg->cigar.empty() && nxt->cigar.empty())
    //            {
    //                break;                       // 不删除任何一条
    //            }

    //            /* 选出较短者作为待删对象 */
    //            SegPtr victim = (seg->length <= nxt->length) ? seg : nxt;
    //            SegPtr survivor = (victim == seg) ? nxt : seg;

    //            /* ---- 真正删除 ---- */
    //            Segment::deleteSegment(victim);  // 已包含 unlink + 从 block 的 anchors 擦除

    //            /* 维护采样表 */
    //            //invalidateSampling(victim->start,
    //                //victim->start + victim->length);

    //            /* 继续比较 survivor 与下一个 */
    //            seg = survivor;
    //            nxt = seg->primary_path.next.load(std::memory_order_relaxed);
    //        }

    //        /* 无重叠，正常前进 */
    //        seg = nxt;
    //    }
    //}




    /* =============================================================
     * 2. Block factories
     * ===========================================================*/
    BlockPtr Block::create(std::size_t hint)
    {
        auto bp = std::make_shared<Block>();
        bp->block_id = g_next_block_id.fetch_add(1, std::memory_order_relaxed);
        bp->anchors.reserve(hint);
        return bp;
    }

    BlockPtr Block::createEmpty(const SpeciesName& ref_species,
                                const ChrName& ref_chr,
                                std::size_t hint)
    {
        if (ref_species.empty() || ref_chr.empty()) {
            throw std::invalid_argument(
                "Block::createEmpty requires a non-empty reference key");
        }
        auto bp = Block::create(hint);
        bp->ref_species = ref_species;
        bp->ref_chr = ref_chr;
        return bp;
    }

    std::pair<SegPtr, SegPtr> Block::createSegmentPair(const Match& match,
        const SpeciesName& ref_name,
        const SpeciesName& qry_name,
        const ChrName& ref_chr,
        const ChrName& qry_chr,
        const BlockPtr& blk)
    {
        // Create ref segment
        uint_t match_len = match.match_len();
        SegPtr ref_seg = Segment::create(match.ref_start, match_len,
            Strand::FORWARD, Cigar_t{ cigarToInt('M', match_len) }, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::create(match.qry_start, match_len,
            match.strand(), Cigar_t{ cigarToInt('M', match_len) }, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // The block is not published yet; no lock is required here.
        blk->ref_species = ref_name;
        blk->anchors.emplace(SpeciesChrPair{ref_name, ref_chr}, ref_seg);
        blk->anchors.emplace(SpeciesChrPair{qry_name, qry_chr}, qry_seg);
        return { ref_seg, qry_seg };
    }

    std::pair<SegPtr, SegPtr> Block::createSegmentPair(const Anchor& anchor,
        const SpeciesName& ref_name,
        const SpeciesName& qry_name,
        const ChrName& ref_chr,
        const ChrName& qry_chr,
        const BlockPtr& blk)
    {
        // Create ref segment
        SegPtr ref_seg = Segment::create(anchor.ref_start, anchor.ref_len,
            Strand::FORWARD, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // Create qry segment
        SegPtr qry_seg = Segment::create(anchor.qry_start, anchor.qry_len,
            anchor.strand, anchor.cigar, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);

        // The block is not published yet; no lock is required here.
        blk->ref_species = ref_name;
        blk->anchors.emplace(SpeciesChrPair{ref_name, ref_chr}, ref_seg);
        blk->anchors.emplace(SpeciesChrPair{qry_name, qry_chr}, qry_seg);
        return { ref_seg, qry_seg };
    }

    std::pair<SegPtr, SegPtr> Block::createSegmentPair(Anchor& anchor,
        const SpeciesName& ref_name,
        const SpeciesName& qry_name,
        const ChrName& ref_chr,
        const ChrName& qry_chr,
        const BlockPtr& blk)
    {
        SegPtr ref_seg = Segment::create(anchor.ref_start, anchor.ref_len,
            Strand::FORWARD, Cigar_t{}, AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);
        SegPtr qry_seg = Segment::create(anchor.qry_start, anchor.qry_len,
            anchor.strand, std::move(anchor.cigar), AlignRole::PRIMARY,
            SegmentRole::SEGMENT, blk);
        blk->ref_species = ref_name;
        blk->anchors.emplace(SpeciesChrPair{ref_name, ref_chr}, ref_seg);
        blk->anchors.emplace(SpeciesChrPair{qry_name, qry_chr}, qry_seg);
        return {ref_seg, qry_seg};
    }
    
    //void Block::removeAllSegments() {
    //    std::unique_lock lk(rw);
    //    
    //    // 从所有链表中解除链接并清理anchors
    //    for (auto& [species_chr, segment] : anchors) {
    //        if (segment) {
    //            Segment::unlinkSegment(segment);
    //            segment->parent_block.reset();
    //            delete segment;
    //        }
    //    }
    //    anchors.clear();
    //}
    void Block::removeAllSegments()
    {
        //std::unique_lock lk(rw);                     // 独占 Block

        std::vector<SegPtr> seg_list;
        seg_list.reserve(anchors.size());

        for (auto& [_, seg] : anchors)       // 只读，绝不修改 anchors
            if (seg) seg_list.push_back(seg);

        anchors.clear();                              // 现在可以一次性清空

        ///* 真正断链 + 释放在容器之外进行，
        //   即使过程中触发再修改 anchors，也不会再有活迭代器。 */
        for (SegPtr seg : seg_list) {
            Segment::unlinkSegment(seg);
            seg->parent_block.reset();
        }

    }

    
    void Block::removeSegmentsBySpecies(const SpeciesName& species) {
        std::unique_lock lk(rw);
        
        std::vector<SpeciesChrPair> to_remove;
        
        // 收集需要删除的entries
        for (const auto& [species_chr, segment] : anchors) {
            if (species_chr.first == species) {
                to_remove.emplace_back(species_chr);
                if (segment) {
                    Segment::unlinkSegment(segment);
                    segment->parent_block.reset();
                }
            }
        }
        
        // 从anchors中移除
        for (const auto& key : to_remove) {
            anchors.erase(key);
        }
    }
    
    size_t Block::removeSegments(const SpeciesName& species, const ChrName& chr) {
        std::unique_lock lk(rw);

        const SpeciesChrPair key{species, chr};
        auto [first, last] = anchors.equal_range(key);
        size_t removed = 0;
        for (auto it = first; it != last; ++it) {
            if (it->second) {
                Segment::unlinkSegment(it->second);
                it->second->parent_block.reset();
            }
            ++removed;
        }
        anchors.erase(first, last);
        return removed;
    }

} // namespace RaMesh
