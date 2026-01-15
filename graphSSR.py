#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SSR Finder - Graph Theory Edition (v7.9 - Final Merging Logic Fix)

This definitive version corrects a fundamental flaw in the SSR reconstruction
logic by implementing a "merge-then-filter" strategy. It ensures that all SSR
fragments, regardless of initial length, are first merged with adjacent,
homologous fragments before the final length filter is applied. This correctly
reconstructs complete SSR loci that have been interrupted by mutations.

Key Scientific Method:
1.  De Bruijn Graph Construction & Cycle Search: Identifies core SSR motifs.
2.  Positional Clustering: Identifies the initial boundaries of all potential
    SSR fragments.
3.  Intelligent Fragment Merging (Corrected): All adjacent fragments sharing
    the same canonical motif are merged into longer candidate SSRs.
4.  Final Length Filtration: The minimum length requirement is applied only AFTER
    the merging process is complete, ensuring that complete, mutation-interrupted
    SSRs are not prematurely discarded.
"""
import time
import os
import argparse
import heapq
import math
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict, Counter

ADAPTIVE_QUANTILE = 0.95

# --- 依赖检查 ---
try:
    import pyfastx
except ImportError:
    pyfastx = None

if pyfastx is not None and os.environ.get("PYFASTX_DISABLE") == "1":
    pyfastx = None


def iter_fasta(path):
    if pyfastx is not None:
        try:
            for record in pyfastx.Fasta(path):
                yield record.name, record.seq
            return
        except Exception:
            pass
    seq_id = None
    chunks = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_id is not None:
                    yield seq_id, ''.join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if seq_id is not None:
            yield seq_id, ''.join(chunks)

# ==============================================================================
#  1. 图构建与分析模块
# ==============================================================================
class DeBruijnGraph:
    """一个简单的德布鲁因图实现。"""
    def __init__(self, sequence, k, min_kmer_count=0, max_cycle_len=None):
        self.k = k
        self.sequence = sequence
        self.min_kmer_count = min_kmer_count
        self.max_cycle_len = max_cycle_len
        self.nodes = Counter()
        self.graph = defaultdict(Counter)
        self.kmer_positions = defaultdict(list)
        self.branch_ratio_threshold = None
        self._build()
        self.branch_ratio_threshold = self._auto_branch_ratio_threshold()

    @staticmethod
    def _median(values):
        if not values:
            return 0.0
        values = sorted(values)
        n = len(values)
        mid = n // 2
        if n % 2 == 1:
            return float(values[mid])
        return (values[mid - 1] + values[mid]) / 2.0

    @staticmethod
    def _quantile(values, q):
        if not values:
            return None
        values = sorted(values)
        n = len(values)
        if n == 1:
            return float(values[0])
        pos = q * (n - 1)
        lower = int(math.floor(pos))
        upper = int(math.ceil(pos))
        if lower == upper:
            return float(values[lower])
        weight = pos - lower
        return values[lower] + (values[upper] - values[lower]) * weight

    def _auto_branch_ratio_threshold(self):
        ratios = []
        for _, successors in self.graph.items():
            if len(successors) < 3:
                continue
            ranked = sorted(successors.values(), reverse=True)
            second_weight = ranked[1]
            third_weight = ranked[2]
            if second_weight <= 0:
                continue
            ratios.append(third_weight / second_weight)
        if not ratios:
            return None
        med = self._median(ratios)
        mad = self._median([abs(r - med) for r in ratios])
        if mad > 0:
            threshold = med + 3 * mad
        else:
            threshold = self._quantile(ratios, ADAPTIVE_QUANTILE)
        if threshold is None:
            return None
        threshold = max(0.0, min(1.0, threshold))
        return threshold

    def _auto_min_kmer_count(self, counts):
        if not counts:
            return 1
        values = list(counts.values())
        med = self._median(values)
        deviations = [abs(v - med) for v in values]
        mad = self._median(deviations)
        if mad > 0:
            threshold = med + 3 * mad
        else:
            threshold = self._quantile(values, ADAPTIVE_QUANTILE)
        if threshold is None:
            threshold = 1
        return max(1, int(math.ceil(threshold)))

    def _auto_cluster_gap(self, position_lists, motif_len, max_merge_dist):
        last_pos = None
        gaps = []
        for pos in heapq.merge(*position_lists):
            if last_pos is not None and pos == last_pos:
                continue
            if last_pos is not None:
                gaps.append(pos - last_pos)
            last_pos = pos
        if not gaps:
            gap = motif_len
        else:
            med = self._median(gaps)
            mad = self._median([abs(g - med) for g in gaps])
            if mad > 0:
                gap = math.ceil(med + 3 * mad)
            else:
                q = self._quantile(gaps, ADAPTIVE_QUANTILE)
                if q is None:
                    gap = motif_len
                else:
                    gap = math.ceil(q)
        gap = max(1, gap)
        if max_merge_dist is not None and max_merge_dist > 0:
            gap = min(gap, max_merge_dist)
        return gap

    def _adaptive_cycle_len(self, start_node):
        positions = self.kmer_positions.get(start_node)
        if not positions or len(positions) < 3:
            return self.k + 7
        gaps = []
        last_pos = None
        for pos in positions:
            if last_pos is not None:
                gap = pos - last_pos
                if gap > 0:
                    gaps.append(gap)
            last_pos = pos
        if len(gaps) < 2:
            return self.k + 7
        med = self._median(gaps)
        mad = self._median([abs(g - med) for g in gaps])
        if mad > 0:
            est = med + 3 * mad
        else:
            est = self._quantile(gaps, ADAPTIVE_QUANTILE)
        if est is None or est <= 0:
            return self.k + 7
        return max(1, int(math.ceil(est)))

    def _build(self):
        """从序列构建德布鲁因图。"""
        if len(self.sequence) < self.k:
            return

        if self.min_kmer_count <= 0:
            counts = Counter()
            for i in range(len(self.sequence) - self.k + 1):
                kmer = self.sequence[i : i + self.k]
                if 'N' in kmer:
                    continue
                counts[kmer] += 1
            self.nodes = counts
            self.min_kmer_count = self._auto_min_kmer_count(counts)
            for i in range(len(self.sequence) - self.k + 1):
                kmer = self.sequence[i : i + self.k]
                if 'N' in kmer:
                    continue
                if self.nodes[kmer] < self.min_kmer_count:
                    continue
                self.kmer_positions[kmer].append(i)

                if i < len(self.sequence) - self.k:
                    next_kmer = self.sequence[i + 1 : i + self.k + 1]
                    if 'N' not in next_kmer and self.nodes[next_kmer] >= self.min_kmer_count:
                        self.graph[kmer][next_kmer] += 1
        elif self.min_kmer_count > 1:
            for i in range(len(self.sequence) - self.k + 1):
                kmer = self.sequence[i : i + self.k]
                if 'N' in kmer:
                    continue
                self.nodes[kmer] += 1

            for i in range(len(self.sequence) - self.k + 1):
                kmer = self.sequence[i : i + self.k]
                if 'N' in kmer:
                    continue
                if self.nodes[kmer] < self.min_kmer_count:
                    continue
                self.kmer_positions[kmer].append(i)

                if i < len(self.sequence) - self.k:
                    next_kmer = self.sequence[i + 1 : i + self.k + 1]
                    if 'N' not in next_kmer and self.nodes[next_kmer] >= self.min_kmer_count:
                        self.graph[kmer][next_kmer] += 1
        else:
            for i in range(len(self.sequence) - self.k + 1):
                kmer = self.sequence[i : i + self.k]
                if 'N' in kmer:
                    continue

                self.nodes[kmer] += 1
                self.kmer_positions[kmer].append(i)

                if i < len(self.sequence) - self.k:
                    next_kmer = self.sequence[i + 1 : i + self.k + 1]
                    if 'N' not in next_kmer:
                        self.graph[kmer][next_kmer] += 1

    def find_ssr_fragments(self, min_avg_coverage, max_merge_dist):
        """在图中寻找代表SSR的高覆盖度循环路径，返回所有碎片。"""
        ssr_fragments = []
        processed_cycles = set()

        sorted_nodes = sorted(self.nodes.items(), key=lambda item: item[1], reverse=True)

        for start_node, coverage in sorted_nodes:
            if coverage < min_avg_coverage:
                break
            
            cycle = self._find_cycle(start_node)
            if not cycle:
                continue

            motif_seq = self._get_motif_from_cycle(cycle)
            motif_len = len(motif_seq)
            if motif_len == 0:
                continue

            canonical_motif = self._get_canonical_motif(motif_seq)
            if canonical_motif in processed_cycles:
                continue
            processed_cycles.add(canonical_motif)

            avg_cov = sum(self.nodes[n] for n in cycle) / len(cycle)
            if avg_cov < min_avg_coverage:
                continue

            position_iters = [self.kmer_positions[node] for node in cycle if node in self.kmer_positions]
            if not position_iters:
                continue

            gap = self._auto_cluster_gap(position_iters, motif_len, max_merge_dist)
            merged_positions = heapq.merge(*position_iters)
            clusters = []
            current_start = None
            current_last = None
            last_pos = None
            for pos in merged_positions:
                if last_pos is not None and pos == last_pos:
                    continue
                if current_start is None:
                    current_start = pos
                    current_last = pos
                elif pos - last_pos > gap:
                    clusters.append((current_start, current_last))
                    current_start = pos
                    current_last = pos
                else:
                    current_last = pos
                last_pos = pos
            if current_start is not None:
                clusters.append((current_start, current_last))

            for start_pos, end_pos in clusters:
                end_pos = end_pos + self.k
                ssr_fragments.append({
                    'Start': start_pos + 1,
                    'End': end_pos,
                    'Inferred_Motif': canonical_motif,
                    'Avg_Coverage': avg_cov
                })
        
        return ssr_fragments

    def _find_cycle(self, start_node):
        """使用深度优先搜索寻找一个简单的环路。"""
        if self.max_cycle_len is None:
            max_len = self._adaptive_cycle_len(start_node)
        else:
            max_len = self.max_cycle_len
        cycle = self._find_cycle_with_branches(start_node, max_len, 2)
        if cycle:
            return cycle
        return self._find_cycle_with_branches(
            start_node,
            max_len,
            2,
            adaptive_ratio=self.branch_ratio_threshold
        )

    def _select_successors(self, successors, max_branches, adaptive_ratio=None):
        ranked = sorted(successors.items(), key=lambda item: (-item[1], item[0]))
        limit = max_branches
        if adaptive_ratio is not None and len(ranked) >= 3:
            second_weight = ranked[1][1]
            third_weight = ranked[2][1]
            if second_weight > 0 and third_weight >= adaptive_ratio * second_weight:
                limit = max(limit, 3)
        return ranked[:limit]

    def _find_cycle_with_branches(self, start_node, max_len, max_branches, adaptive_ratio=None):
        stack = [(start_node, [start_node], {start_node})]
        while stack:
            current_node, path, visited_on_path = stack.pop()
            if len(path) > max_len:
                continue
            successors = self.graph.get(current_node)
            if not successors:
                continue
            ranked = self._select_successors(successors, max_branches, adaptive_ratio=adaptive_ratio)
            for next_node, _ in ranked:
                if next_node in visited_on_path:
                    try:
                        cycle_start_index = path.index(next_node)
                        return path[cycle_start_index:]
                    except ValueError:
                        continue
                if len(path) < max_len:
                    new_path = path + [next_node]
                    new_visited = set(visited_on_path)
                    new_visited.add(next_node)
                    stack.append((next_node, new_path, new_visited))
        return None

    def _get_motif_from_cycle(self, cycle):
        """从环路路径重建基元序列。"""
        if not cycle: return ""
        return "".join(node[0] for node in cycle)

    def _get_canonical_motif(self, motif):
        motif_rc = str.maketrans('ATCG', 'TAGC')
        rc_motif = motif.translate(motif_rc)[::-1]
        candidates = []
        for m in [motif, rc_motif]:
            if not m: continue
            for i in range(len(m)):
                candidates.append(m[i:] + m[:i])
        return min(candidates) if candidates else ""

# ==============================================================================
#  4. 主分析流程
# ==============================================================================
def find_ssrs_in_sequence(sequence, sequence_id, k, min_len, min_cov, max_merge_dist,
                          min_kmer_count=0, max_cycle_len=None, verbose=True):
    """分析单条序列的主函数。"""
    if verbose:
        print(f"\n--- 正在处理序列: {sequence_id} (长度: {len(sequence)} bp) ---")
        print(f"  - 步骤 1: 正在构建 k={k} 的德布鲁因图...")
    graph = DeBruijnGraph(
        sequence,
        k,
        min_kmer_count=min_kmer_count,
        max_cycle_len=max_cycle_len
    )
    
    if verbose:
        print(f"  - 步骤 2: 正在图中寻找所有SSR碎片...")
    # 先不过滤长度，找到所有可能的碎片
    ssr_fragments = graph.find_ssr_fragments(min_avg_coverage=min_cov, max_merge_dist=max_merge_dist)
    
    if verbose:
        print(f"  - 步骤 3: 正在合并和修剪 {len(ssr_fragments)} 个候选碎片...")
    if not ssr_fragments:
        if verbose:
            print(f"  - 步骤 4: 未找到任何最终SSR。")
        return []

    ssr_fragments.sort(key=lambda x: x['Start'])

    merged_ssrs = []
    if ssr_fragments:
        current_ssr = ssr_fragments[0]
        for i in range(1, len(ssr_fragments)):
            next_ssr = ssr_fragments[i]
            dist = next_ssr['Start'] - current_ssr['End']

            if next_ssr['Inferred_Motif'] == current_ssr['Inferred_Motif'] and dist <= max_merge_dist:
                current_ssr['End'] = max(current_ssr['End'], next_ssr['End'])
            else:
                merged_ssrs.append(current_ssr)
                current_ssr = next_ssr
        merged_ssrs.append(current_ssr)

    final_ssrs = []
    for ssr in merged_ssrs:
        ssr['Size'] = ssr['End'] - ssr['Start'] + 1
        
        # 只有在合并之后，才进行长度过滤
        if ssr['Size'] >= min_len:
            motif_len = len(ssr['Inferred_Motif'])
            if motif_len == 0: continue
            
            # 精确修剪
            new_size = (ssr['Size'] // motif_len) * motif_len
            if new_size < motif_len:
                continue
            ssr['Size'] = new_size
            ssr['End'] = ssr['Start'] + new_size - 1
            ssr['SSR_Sequence'] = sequence[ssr['Start'] - 1 : ssr['End']]
            ssr['sequence_id'] = sequence_id
            final_ssrs.append(ssr)
            
    if verbose:
        print(f"  - 步骤 4: 已找到 {len(final_ssrs)} 个最终SSR。")
    return final_ssrs

# ==============================================================================
#  5. 结果写入与主执行模块
# ==============================================================================
def write_results_to_file(all_ssrs, output_filename):
    if not all_ssrs:
        print("\n所有序列中均未找到SSR，不生成输出文件。")
        return
    print(f"\n--- 正在写入总文件 ---")
    
    all_ssrs.sort(key=lambda x: (x['sequence_id'], x['Start']))
    for i, ssr in enumerate(all_ssrs):
        ssr['SSR_nr'] = i + 1

    with open(output_filename, 'w') as f:
        header = "Sequence_ID\tSSR_nr\tInferred_Motif\tSSR_Sequence\tSize\tStart\tEnd\tAvg_Coverage\n"
        f.write(header)
        for ssr in all_ssrs:
            line = (
                f"{ssr['sequence_id']}\t{ssr['SSR_nr']}\t{ssr['Inferred_Motif']}\t"
                f"{ssr['SSR_Sequence']}\t{ssr['Size']}\t{ssr['Start']}\t"
                f"{ssr['End']}\t{ssr['Avg_Coverage']:.2f}"
            )
            f.write(line + "\n")
    
    print(f"  - 详细结果已保存到: {output_filename}")


def _process_sequence(args):
    seq_id, seq, k, min_len, min_cov, max_merge_dist = args
    return find_ssrs_in_sequence(
        seq,
        seq_id,
        k,
        min_len,
        min_cov,
        max_merge_dist,
        min_kmer_count=0,
        max_cycle_len=None,
        verbose=False
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SSR Finder - Graph Theory Edition: 基于德布鲁因图的SSR识别工具。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_fasta", help="输入的FASTA格式的基因组序列文件。")
    parser.add_argument("output_prefix", help="输出文件的前缀 (例如: 'my_genome_ssrs')。")
    parser.add_argument("-k", "--kmer_size", type=int, default=21,
                        help="用于构建德布鲁因图的k-mer大小 (默认: 21)。\n"
                             "这是最关键的参数，决定了图的分辨率。")
    parser.add_argument("--min_len", type=int, default=30,
                        help="报告SSR的最小长度 (默认: 30 bp)。")
    parser.add_argument("--min_cov", type=float, default=2.0,
                        help="报告SSR环路的最小平均k-mer覆盖度 (默认: 2.0)。\n"
                             "用于将重复区域与单拷贝区域区分开。")
    parser.add_argument("--max_merge_dist", type=int, default=25,
                        help="合并邻近SSR碎片的最大允许间隔距离 (默认: 25 bp)。")
    parser.add_argument("--threads", type=int, default=1,
                        help="并行处理的进程数 (默认: 1；0 表示自动使用CPU核数)。")
    parser.add_argument("--quiet", action='store_true',
                        help="关闭逐序列日志输出。")
    args = parser.parse_args()

    try:
        start_time = time.time()
        total_sequences, all_ssrs_found = 0, []
        
        print("\n" + "="*50)
        print(f"开始分析文件: {args.input_fasta}")
        print("="*50)
        
        if pyfastx is None:
            print("警告: 未检测到 pyfastx，使用内置FASTA解析器。")

        threads = args.threads
        if threads == 0:
            threads = os.cpu_count() or 1
        if threads < 1:
            print("错误: --threads 必须 >= 1 或 0(自动)。")
            exit(1)

        if threads == 1:
            for seq_id, seq in iter_fasta(args.input_fasta):
                total_sequences += 1
                processed_ssrs = find_ssrs_in_sequence(
                    seq, seq_id,
                    args.kmer_size, args.min_len, args.min_cov, args.max_merge_dist,
                    min_kmer_count=0,
                    max_cycle_len=None,
                    verbose=not args.quiet
                )
                all_ssrs_found.extend(processed_ssrs)
        else:
            sequences = list(iter_fasta(args.input_fasta))
            total_sequences = len(sequences)
            print(f"并行处理 {total_sequences} 条序列 (threads={threads})...")
            task_args = [
                (seq_id, seq, args.kmer_size, args.min_len, args.min_cov, args.max_merge_dist)
                for seq_id, seq in sequences
            ]
            with ProcessPoolExecutor(max_workers=threads) as executor:
                for processed_ssrs in executor.map(_process_sequence, task_args):
                    all_ssrs_found.extend(processed_ssrs)

        output_filename = f"{args.output_prefix}.ssrs"
        write_results_to_file(all_ssrs_found, output_filename)
        end_time = time.time()
        
        print("\n" + "="*50)
        print("--- 分析全部完成 ---")
        print(f"共处理了 {total_sequences} 条序列，总计发现 {len(all_ssrs_found)} 个SSR。")
        print(f"结果已保存到 {output_filename}")
        print(f"总耗时: {end_time - start_time:.2f} 秒")
        print("="*50)

    except FileNotFoundError:
        print(f"错误: 输入文件未找到 -> {args.input_fasta}")
        exit(1)
    except Exception as e:
        print(f"发生了一个意外错误: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
