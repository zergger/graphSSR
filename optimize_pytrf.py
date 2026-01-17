#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import argparse
import random
import textwrap
import shutil
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    print("错误: 本脚本需要 'numpy' 库来进行统计计算。")
    print("请运行: pip install numpy")
    exit(1)

try:
    from scipy import stats
except ImportError:
    print("错误: 本脚本需要 'scipy' 库来进行统计检验。")
    print("请运行: pip install scipy")
    exit(1)


# --- 模块1: SSR进化模拟器 (与TRF版本相同，无需改动) ---
def generate_combined_dataset(motifs_to_test, repeats_map, flank_len, snp_rate, indel_rate, work_dir, replicate_id=0):
    fasta_path = os.path.join(work_dir, f"combined_rep_{replicate_id}.fasta")
    truth_path = os.path.join(work_dir, f"combined_truth_rep_{replicate_id}.tsv")
    bases = ['A', 'T', 'C', 'G']

    with open(fasta_path, 'w') as fasta_f, open(truth_path, 'w') as truth_f:
        truth_f.write("SequenceID\tMotifType\tTrue_Start\tTrue_End\n")
        
        for motif_type, motif in motifs_to_test.items():
            repeats = repeats_map[motif_type]
            left_flank_perfect = "".join(random.choices(bases, k=flank_len))
            right_flank_perfect = "".join(random.choices(bases, k=flank_len))
            ssr_core_perfect = motif * repeats
            left_flank_mut, core_mut, right_flank_mut = [], [], []

            def evolve(sequence, target_list):
                for base in sequence:
                    rand_val = random.random()
                    if rand_val < indel_rate:
                        if random.choice(['ins', 'del']) == 'ins':
                            indel_len = random.randint(1, 2)
                            insertion = random.choices(bases, k=indel_len)
                            target_list.extend(insertion)
                            target_list.append(base)
                    elif rand_val < indel_rate + snp_rate:
                        target_list.append(random.choice([b for b in bases if b != base]))
                    else:
                        target_list.append(base)

            evolve(left_flank_perfect, left_flank_mut)
            evolve(ssr_core_perfect, core_mut)
            evolve(right_flank_perfect, right_flank_mut)

            evolved_seq_str = "".join(left_flank_mut + core_mut + right_flank_mut)
            new_ssr_start = len(left_flank_mut) + 1
            new_ssr_end = len(left_flank_mut) + len(core_mut)
            
            seq_id = f"sim_{motif_type}_rep{replicate_id}"
            fasta_f.write(f">{seq_id}\n" + "\n".join(textwrap.wrap(evolved_seq_str, 60)) + "\n")
            truth_f.write(f"{seq_id}\t{motif_type}\t{new_ssr_start}\t{new_ssr_end}\n")
            
    return fasta_path, truth_path


def build_sim_metadata(motifs_to_test, repeats_map, flank_len, snp_rate, indel_rate, replicate_id, seed):
    params = {
        'motifs_to_test': motifs_to_test,
        'repeats_map': repeats_map,
        'flank_len': flank_len,
        'snp_rate': snp_rate,
        'indel_rate': indel_rate,
        'replicate_id': replicate_id,
        'seed': seed,
    }
    signature = json.dumps(params, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return {'signature': signature, 'params': params}


def load_sim_metadata(meta_path):
    try:
        with open(meta_path, 'r') as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def write_sim_metadata(meta_path, meta):
    tmp_path = f"{meta_path}.tmp"
    with open(tmp_path, 'w') as f:
        json.dump(meta, f, ensure_ascii=True, sort_keys=True, indent=2)
    os.replace(tmp_path, meta_path)


def read_fasta_lengths(fasta_path):
    lengths = {}
    seq_id = None
    current_len = 0
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if seq_id is not None:
                        lengths[seq_id] = current_len
                    seq_id = line[1:].strip().split()[0]
                    current_len = 0
                else:
                    current_len += len(line.strip())
        if seq_id is not None:
            lengths[seq_id] = current_len
    except OSError:
        return {}
    return lengths


def parse_pytrf_predictions(pytrf_output_file):
    predictions_by_seq = defaultdict(list)
    if not os.path.exists(pytrf_output_file) or os.path.getsize(pytrf_output_file) == 0:
        return predictions_by_seq
    try:
        pred_df = pd.read_csv(pytrf_output_file, sep='\t', header=None)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        pred_df = None

    if pred_df is not None and not pred_df.empty:
        for row in pred_df.itertuples(index=False, name=None):
            if len(row) < 3:
                continue
            seq_id = row[0]
            try:
                start = int(row[1]) + 1  # pytrf start is 0-based
                end = int(row[2])
            except (ValueError, TypeError):
                continue
            period = None
            if len(row) > 4:
                try:
                    period = int(float(row[4]))
                except (ValueError, TypeError):
                    period = None
            pred = {'start': start, 'end': end, 'period': period}
            predictions_by_seq[seq_id].append(pred)
    return predictions_by_seq


def compute_counts_by_motif(predictions_by_seq, truth_data, motifs_to_test, match_mode, overlap_threshold):
    counts_by_motif = {mtype: {'tp': 0, 'fp': 0, 'fn': 0} for mtype in motifs_to_test}

    truths_by_motif = defaultdict(list)
    for t in truth_data:
        truths_by_motif[t['MotifType']].append(t)

    for motif_type in motifs_to_test:
        tp = fp = fn = 0
        truths_for_motif = truths_by_motif.get(motif_type, [])
        truth_by_seq = {t['SequenceID']: t for t in truths_for_motif}
        truth_seq_ids = set(truth_by_seq.keys())
        pred_seq_ids_for_motif = {sid for sid in predictions_by_seq if motif_type in sid}
        all_seq_ids = truth_seq_ids.union(pred_seq_ids_for_motif)
        required_period = len(motifs_to_test[motif_type]) if match_mode == "interval+period" else None

        for seq_id in all_seq_ids:
            preds = predictions_by_seq.get(seq_id, [])
            truth = truth_by_seq.get(seq_id)

            if truth:
                is_matched = False
                true_start, true_end = truth['True_Start'], truth['True_End']
                for pred in preds:
                    if required_period is not None and pred.get('period') != required_period:
                        continue
                    pred_start, pred_end = pred['start'], pred['end']
                    intersection = max(0, min(true_end, pred_end) - max(true_start, pred_start) + 1)
                    union = (true_end - true_start + 1) + (pred_end - pred_start + 1) - intersection
                    overlap = intersection / union if union > 0 else 0
                    if overlap >= overlap_threshold:
                        is_matched = True
                        break
                if is_matched:
                    tp += 1
                    fp += len(preds) - 1
                else:
                    fn += 1
                    fp += len(preds)
            else:
                fp += len(preds)

        counts_by_motif[motif_type] = {'tp': tp, 'fp': fp, 'fn': fn}

    return counts_by_motif


def compute_metrics_from_counts(counts_by_motif, motifs_to_test):
    f1_scores = {}
    precisions = []
    recalls = []
    f1s = []
    for motif_type in motifs_to_test:
        counts = counts_by_motif.get(motif_type, {'tp': 0, 'fp': 0, 'fn': 0})
        tp = counts['tp']
        fp = counts['fp']
        fn = counts['fn']
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        f1_scores[motif_type] = f1
        precisions.append(precision)
        recalls.append(recall)
        f1s.append(f1)

    macro_precision = float(np.mean(precisions)) if precisions else 0.0
    macro_recall = float(np.mean(recalls)) if recalls else 0.0
    macro_f1 = float(np.mean(f1s)) if f1s else 0.0
    return f1_scores, macro_precision, macro_recall, macro_f1


def compute_null_overlaps(predictions_by_seq, truth_by_seq, seq_lengths, rng):
    overlaps = []
    for seq_id, preds in predictions_by_seq.items():
        truth = truth_by_seq.get(seq_id)
        seq_len = seq_lengths.get(seq_id)
        if not truth or not seq_len:
            continue
        true_start, true_end = truth['True_Start'], truth['True_End']
        for pred in preds:
            length = pred['end'] - pred['start'] + 1
            if length <= 0 or length > seq_len:
                continue
            rand_start = rng.randint(1, seq_len - length + 1)
            rand_end = rand_start + length - 1
            intersection = max(0, min(true_end, rand_end) - max(true_start, rand_start) + 1)
            union = (true_end - true_start + 1) + length - intersection
            if union > 0:
                overlaps.append(intersection / union)
    return overlaps


# --- 模块2: 评估器 (适配pytrf输出) ---
def evaluate_f1_by_motif(pytrf_output_file, truth_file, motifs_to_test, match_mode="interval", overlap_threshold=0.5):
    f1_scores = {mtype: 0.0 for mtype in motifs_to_test}
    
    if not os.path.exists(pytrf_output_file) or os.path.getsize(pytrf_output_file) == 0:
        return f1_scores
        
    truth_data = pd.read_csv(truth_file, sep='\t').to_dict('records')

    predictions_by_seq = parse_pytrf_predictions(pytrf_output_file)
    counts_by_motif = compute_counts_by_motif(
        predictions_by_seq,
        truth_data,
        motifs_to_test,
        match_mode,
        overlap_threshold
    )
    f1_scores, _, _, _ = compute_metrics_from_counts(counts_by_motif, motifs_to_test)
    return f1_scores

# --- 模块3: 主控流程 ---
def main():
    parser = argparse.ArgumentParser(description="高级自动化pytrf参数一体化网格搜索与统计优化流程 (重复模拟法)。", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--num_replicates", type=int, default=30, help="为进行统计检验而执行的独立模拟重复次数。")
    # pytrf特有的参数
    parser.add_argument("--min_identity_range", default="0.7,0.8,0.9", help="测试的最小相似度(min-identity)值范围, 以逗号分隔。")
    parser.add_argument("--max_error_range", default="2,3,4", help="测试的最大连续错误数(max-continuous-error)值范围, 以逗号分隔。")
    parser.add_argument("--min_seed_repeat_range", default="3,4,5", help="测试的最小种子重复数(min-seed-repeat)值范围, 以逗号分隔。")
    # 模拟参数
    parser.add_argument("--snp_rate", type=float, default=0.02, help="模拟的SNP频率。")
    parser.add_argument("--indel_rate", type=float, default=0.01, help="模拟的Indel频率。")
    parser.add_argument("--seed", type=int, default=None, help="随机种子基值(可复现); 如设置，则每次重复使用 seed + replicate_id。")
    parser.add_argument("--output_csv", default="pytrf_grid_search_results.csv", help="保存网格搜索结果的CSV文件名。")
    parser.add_argument("--keep_files", action='store_true', help="保留所有中间生成的模拟和pytrf结果文件。")
    parser.add_argument("--reuse_sim", action='store_true', help="复用已存在的模拟数据(忽略参数匹配检查)。")
    parser.add_argument("--force_regen", action='store_true', help="强制重新生成模拟数据。")
    parser.add_argument("--match_mode", choices=["interval", "interval+period"], default="interval", help="评估模式: interval=仅区间重叠; interval+period=区间重叠且period匹配。")
    parser.add_argument("--max_workers", type=int, default=1, help="并行运行pytrf的最大并发数(每个replicate内)。")
    parser.add_argument("--overlap_threshold", type=float, default=0.5, help="区间重叠阈值(IoU)，用于未启用阈值校准时。")
    parser.add_argument("--calibrate_threshold", action='store_true', help="使用校准集+零假设分布自动选择重叠阈值。")
    parser.add_argument("--calibration_fraction", type=float, default=0.2, help="校准集占比(0-1)，仅在未指定 --calibration_replicates 时生效。")
    parser.add_argument("--calibration_replicates", type=int, default=None, help="用于校准阈值的重复次数。")
    parser.add_argument("--threshold_grid", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9", help="阈值扫描列表(逗号分隔)。")
    parser.add_argument("--target_precision", type=float, default=0.9, help="阈值选择的精确率下限(0-1)。")
    parser.add_argument("--null_quantile", type=float, default=0.95, help="零假设重叠分布分位数(0-1)，用于设定t_null。")
    args = parser.parse_args()

    if args.reuse_sim and args.force_regen:
        print("错误: --reuse_sim 与 --force_regen 不能同时使用。"); return
    if args.max_workers < 1:
        print("错误: --max_workers 必须 >= 1。"); return
    if not (0 < args.overlap_threshold <= 1):
        print("错误: --overlap_threshold 必须在 (0, 1] 范围内。"); return
    if not (0 < args.calibration_fraction < 1):
        print("错误: --calibration_fraction 必须在 (0, 1) 范围内。"); return
    if not (0 <= args.target_precision <= 1):
        print("错误: --target_precision 必须在 [0, 1] 范围内。"); return
    if not (0 < args.null_quantile < 1):
        print("错误: --null_quantile 必须在 (0, 1) 范围内。"); return

    identity_values = [float(p) for p in args.min_identity_range.split(',')]
    error_values = [int(p) for p in args.max_error_range.split(',')]
    seed_repeat_values = [int(s) for s in args.min_seed_repeat_range.split(',')]

    motifs_to_test = {'mono': 'A', 'di': 'AT', 'tri': 'AGC', 'tetra': 'AATC', 'penta': 'AGCTC', 'hexa': 'AATGCC'}
    repeats_map = {'mono': 50, 'di': 25, 'tri': 17, 'tetra': 12, 'penta': 10, 'hexa': 8}
    flank_len = 200
    work_dir = "pytrf_grid_search_workspace"
    os.makedirs(work_dir, exist_ok=True)
    
    param_combos = [
        (identity, error, seed_repeat)
        for identity in identity_values
        for error in error_values
        for seed_repeat in seed_repeat_values
    ]

    print("\n" + "="*50)
    print(f"--- 开始pytrf一体化网格搜索 (重复模拟法) ---")
    print(f"总计 {args.num_replicates} 次重复, {len(param_combos)} 组参数。")
    print("="*50)

    def prepare_sim_data(replicate_id):
        fasta_file = os.path.join(work_dir, f"combined_rep_{replicate_id}.fasta")
        truth_file = os.path.join(work_dir, f"combined_truth_rep_{replicate_id}.tsv")
        meta_file = os.path.join(work_dir, f"combined_rep_{replicate_id}.meta.json")

        replicate_seed = args.seed + replicate_id if args.seed is not None else None
        sim_meta = build_sim_metadata(
            motifs_to_test, repeats_map, flank_len,
            args.snp_rate, args.indel_rate, replicate_id, replicate_seed
        )

        reuse_sim = False
        if args.force_regen:
            reuse_sim = False
        elif args.reuse_sim:
            reuse_sim = os.path.exists(fasta_file) and os.path.exists(truth_file)
        else:
            if os.path.exists(fasta_file) and os.path.exists(truth_file) and os.path.exists(meta_file):
                existing_meta = load_sim_metadata(meta_file)
                if existing_meta and existing_meta.get('signature') == sim_meta['signature']:
                    reuse_sim = True

        if reuse_sim:
            print("  - 发现匹配的组合模拟数据集，将重新使用。")
        else:
            print("  - 正在生成组合模拟数据集...")
            if replicate_seed is not None:
                random.seed(replicate_seed)
                np.random.seed(replicate_seed)
            generate_combined_dataset(
                motifs_to_test, repeats_map, flank_len,
                args.snp_rate, args.indel_rate, work_dir, replicate_id=replicate_id
            )
            write_sim_metadata(meta_file, sim_meta)

        return fasta_file, truth_file

    def parse_threshold_grid(grid_str):
        values = []
        for part in grid_str.split(','):
            part = part.strip()
            if not part:
                continue
            values.append(float(part))
        values = sorted(set(values))
        return values

    def build_threshold_output_paths(output_csv):
        base, ext = os.path.splitext(output_csv)
        if not ext:
            base = output_csv
            ext = ".csv"
        curve_path = f"{base}_threshold_curve{ext}"
        meta_path = f"{base}_threshold_meta.json"
        return curve_path, meta_path

    threshold_grid = parse_threshold_grid(args.threshold_grid)
    if not threshold_grid:
        print("错误: --threshold_grid 解析为空。"); return
    if any(t <= 0 or t > 1 for t in threshold_grid):
        print("错误: --threshold_grid 必须在 (0, 1] 范围内。"); return

    if args.calibrate_threshold:
        if args.calibration_replicates is not None:
            calibration_count = args.calibration_replicates
        else:
            calibration_count = max(1, int(args.num_replicates * args.calibration_fraction))
        if calibration_count >= args.num_replicates:
            print("错误: 校准集数量必须小于总重复次数。"); return
        calibration_indices = list(range(calibration_count))
        evaluation_indices = list(range(calibration_count, args.num_replicates))
        print(f"校准集: {len(calibration_indices)} 次重复; 评估集: {len(evaluation_indices)} 次重复。")
    else:
        calibration_indices = []
        evaluation_indices = list(range(args.num_replicates))

    results_dist = defaultdict(list)

    selected_threshold = args.overlap_threshold
    t_null = None

    if args.calibrate_threshold:
        threshold_stats = {
            t: {m: {'tp': 0, 'fp': 0, 'fn': 0} for m in motifs_to_test}
            for t in threshold_grid
        }
        null_overlaps = []
        null_seed_base = args.seed if args.seed is not None else 1

        print("\n--- 阈值校准阶段 ---")
        for i in calibration_indices:
            print(f"\n--- Calibration Replicate {i + 1}/{args.num_replicates} ---")
            fasta_file, truth_file = prepare_sim_data(i)
            truth_data = pd.read_csv(truth_file, sep='\t').to_dict('records')
            truth_by_seq = {t['SequenceID']: t for t in truth_data}
            seq_lengths = read_fasta_lengths(fasta_file)

            print("  - 正在评估所有参数组合(用于阈值校准)...")

            def run_pytrf_for_calibration(j, params):
                identity, error, seed_repeat = params
                pytrf_output_file = os.path.join(work_dir, f"pytrf_combo_{j}_rep_{i}.tsv")

                if not os.path.exists(pytrf_output_file):
                    pytrf_command = [
                        'pytrf', 'findatr', fasta_file,
                        '-o', pytrf_output_file,
                        '-f', 'tsv',
                        '-p', str(identity),
                        '-e', str(error),
                        '-r', str(seed_repeat)
                    ]
                    result = subprocess.run(pytrf_command, capture_output=True, text=True)
                    if result.returncode != 0:
                        stderr_text = (result.stderr or "").strip()
                        if len(stderr_text) > 500:
                            stderr_text = stderr_text[:500] + "..."
                        raise RuntimeError(
                            f"pytrf运行失败: params={(identity, error, seed_repeat)}, returncode={result.returncode}, stderr={stderr_text}"
                        )
                    if not os.path.exists(pytrf_output_file):
                        raise RuntimeError(f"pytrf输出文件未生成: {pytrf_output_file}")

                predictions_by_seq = parse_pytrf_predictions(pytrf_output_file)
                counts_by_threshold = {}
                for t in threshold_grid:
                    counts_by_threshold[t] = compute_counts_by_motif(
                        predictions_by_seq,
                        truth_data,
                        motifs_to_test,
                        args.match_mode,
                        t
                    )
                rng = random.Random(null_seed_base + i * 100000 + j)
                null_ov = compute_null_overlaps(predictions_by_seq, truth_by_seq, seq_lengths, rng)
                return counts_by_threshold, null_ov

            if args.max_workers == 1:
                for j, params in enumerate(param_combos):
                    print(f"\r    - Testing combo {j+1}/{len(param_combos)}...", end="")
                    try:
                        counts_by_threshold, null_ov = run_pytrf_for_calibration(j, params)
                    except Exception as e:
                        print(f"\n!!! 错误: {e}")
                        return
                    for t, counts in counts_by_threshold.items():
                        for motif_type, c in counts.items():
                            total = threshold_stats[t][motif_type]
                            total['tp'] += c['tp']
                            total['fp'] += c['fp']
                            total['fn'] += c['fn']
                    null_overlaps.extend(null_ov)
            else:
                completed = 0
                with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                    futures = {executor.submit(run_pytrf_for_calibration, j, params): params for j, params in enumerate(param_combos)}
                    for future in as_completed(futures):
                        try:
                            counts_by_threshold, null_ov = future.result()
                        except Exception as e:
                            print(f"\n!!! 错误: {e}")
                            return
                        completed += 1
                        print(f"\r    - Testing combo {completed}/{len(param_combos)}...", end="")
                        for t, counts in counts_by_threshold.items():
                            for motif_type, c in counts.items():
                                total = threshold_stats[t][motif_type]
                                total['tp'] += c['tp']
                                total['fp'] += c['fp']
                                total['fn'] += c['fn']
                        null_overlaps.extend(null_ov)
            print("\r    - 完成所有参数组合的评估。           ")

        if null_overlaps:
            t_null = float(np.quantile(null_overlaps, args.null_quantile))
        else:
            t_null = 0.0
            print("警告: 零假设重叠分布为空，将t_null设为0。")

        threshold_records = []
        for t in threshold_grid:
            _, macro_precision, macro_recall, macro_f1 = compute_metrics_from_counts(
                threshold_stats[t], motifs_to_test
            )
            threshold_records.append({
                'threshold': t,
                'macro_precision': macro_precision,
                'macro_recall': macro_recall,
                'macro_f1': macro_f1
            })

        eligible = [r for r in threshold_records if r['threshold'] >= t_null]
        if not eligible:
            eligible = threshold_records[:]
        if args.target_precision > 0:
            eligible_prec = [r for r in eligible if r['macro_precision'] >= args.target_precision]
            candidates = eligible_prec if eligible_prec else eligible
        else:
            candidates = eligible

        candidates = sorted(
            candidates,
            key=lambda r: (r['macro_f1'], r['macro_precision'], r['threshold']),
            reverse=True
        )
        selected_threshold = candidates[0]['threshold']

        curve_path, meta_path = build_threshold_output_paths(args.output_csv)
        pd.DataFrame(threshold_records).to_csv(curve_path, index=False, float_format="%.4f")
        meta = {
            't_null': t_null,
            'selected_threshold': selected_threshold,
            'threshold_grid': threshold_grid,
            'target_precision': args.target_precision,
            'null_quantile': args.null_quantile,
            'calibration_replicates': len(calibration_indices),
            'evaluation_replicates': len(evaluation_indices),
            'match_mode': args.match_mode
        }
        with open(meta_path, 'w') as f:
            json.dump(meta, f, ensure_ascii=True, indent=2, sort_keys=True)

        print("\n--- 阈值校准结果 ---")
        print(f"  - t_null (quantile={args.null_quantile}): {t_null:.4f}")
        print(f"  - 选择的阈值: {selected_threshold:.4f}")
        print(f"  - 阈值曲线已保存: {curve_path}")
        print(f"  - 阈值元数据已保存: {meta_path}")

    print("\n--- 评估阶段 ---")
    for idx, i in enumerate(evaluation_indices, start=1):
        print(f"\n--- Evaluation Replicate {idx}/{len(evaluation_indices)} (原始编号 {i + 1}) ---")
        fasta_file, truth_file = prepare_sim_data(i)

        print("  - 正在评估所有参数组合...")

        def run_pytrf_and_eval(j, params):
            identity, error, seed_repeat = params
            param_tuple = (identity, error, seed_repeat)
            pytrf_output_file = os.path.join(work_dir, f"pytrf_combo_{j}_rep_{i}.tsv")

            if not os.path.exists(pytrf_output_file):
                pytrf_command = [
                    'pytrf', 'findatr', fasta_file,
                    '-o', pytrf_output_file,
                    '-f', 'tsv',
                    '-p', str(identity),
                    '-e', str(error),
                    '-r', str(seed_repeat)
                ]
                result = subprocess.run(pytrf_command, capture_output=True, text=True)
                if result.returncode != 0:
                    stderr_text = (result.stderr or "").strip()
                    if len(stderr_text) > 500:
                        stderr_text = stderr_text[:500] + "..."
                    raise RuntimeError(
                        f"pytrf运行失败: params={param_tuple}, returncode={result.returncode}, stderr={stderr_text}"
                    )
                if not os.path.exists(pytrf_output_file):
                    raise RuntimeError(f"pytrf输出文件未生成: {pytrf_output_file}")

            f1_scores = evaluate_f1_by_motif(
                pytrf_output_file,
                truth_file,
                motifs_to_test,
                match_mode=args.match_mode,
                overlap_threshold=selected_threshold
            )
            return param_tuple, f1_scores

        if args.max_workers == 1:
            for j, params in enumerate(param_combos):
                print(f"\r    - Testing combo {j+1}/{len(param_combos)}...", end="")
                try:
                    param_tuple, f1_scores = run_pytrf_and_eval(j, params)
                except Exception as e:
                    print(f"\n!!! 错误: {e}")
                    return
                results_dist[param_tuple].append(f1_scores)
        else:
            completed = 0
            with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                futures = {executor.submit(run_pytrf_and_eval, j, params): params for j, params in enumerate(param_combos)}
                for future in as_completed(futures):
                    try:
                        param_tuple, f1_scores = future.result()
                    except Exception as e:
                        print(f"\n!!! 错误: {e}")
                        return
                    completed += 1
                    print(f"\r    - Testing combo {completed}/{len(param_combos)}...", end="")
                    results_dist[param_tuple].append(f1_scores)
        print("\r    - 完成所有参数组合的评估。           ")

    # --- 结果分析 ---
    all_results = []
    for param_tuple, f1_scores_per_replicate_list in results_dist.items():
        identity, error, seed_repeat = param_tuple
        a_mean_list = [np.mean(list(f1_dict.values())) for f1_dict in f1_scores_per_replicate_list]
        avg_a_mean = np.mean(a_mean_list)
        ci_lower = np.percentile(a_mean_list, 2.5)
        ci_upper = np.percentile(a_mean_list, 97.5)
        ci_str = f"[{ci_lower:.4f}, {ci_upper:.4f}]"

        current_result = {
            'min_identity': identity, 'max_error': error, 'min_seed_repeat': seed_repeat,
            'f1_mean_avg': avg_a_mean,
            'f1_mean_ci': ci_str,
            'a_mean_dist': a_mean_list,
            'f1_scores_per_replicate': f1_scores_per_replicate_list
        }

        f1_scores_by_motif = defaultdict(list)
        for replicate_scores in f1_scores_per_replicate_list:
            for motif_type, f1_score in replicate_scores.items():
                f1_scores_by_motif[motif_type].append(f1_score)
        
        for motif_type in motifs_to_test.keys():
            f1_list = f1_scores_by_motif[motif_type]
            current_result[f'f1_{motif_type}_avg'] = np.mean(f1_list) if f1_list else 0.0
        
        all_results.append(current_result)
        
    if not all_results:
        print("\n!!! 错误: 未能生成任何结果。"); return

    results_df = pd.DataFrame(all_results)
    
    print("\n" + "="*80)
    print("---                                 pytrf网格搜索完整结果总结                                 ---")
    print("="*80)
    motif_f1_cols = [f'f1_{mtype}_avg' for mtype in motifs_to_test.keys()]
    display_cols = ['min_identity', 'max_error', 'min_seed_repeat'] + motif_f1_cols + ['f1_mean_avg', 'f1_mean_ci']
    print(results_df[display_cols].sort_values(by='f1_mean_avg', ascending=False).to_string(index=False, float_format="%.4f"))

    # --- 最终决策 ---
    if results_df.empty or results_df['f1_mean_avg'].sum() == 0:
        print("\n!!! 警告: 所有参数组合的算术平均F1分数均为0。未能找到最佳参数。")
    else:
        results_df_sorted = results_df.sort_values(by='f1_mean_avg', ascending=False).reset_index(drop=True)
        best_row = results_df_sorted.iloc[0]
        best_params_dict = {
            'min_identity': best_row['min_identity'],
            'max_error': int(best_row['max_error']),
            'min_seed_repeat': int(best_row['min_seed_repeat'])
        }
        
        best_a_mean_dist = best_row['a_mean_dist']
        best_f1_scores_per_replicate = best_row['f1_scores_per_replicate']

        print("\n" + "="*80)
        print("\n--- 最终结论: 推荐使用的最佳统一pytrf参数组合 ---")
        print("  - 最佳参数组合:")
        for key, val in best_params_dict.items():
            print(f"    - {key.replace('_', '-')}: {val}")
            
        mean_a_f1 = np.mean(best_a_mean_dist)
        ci_lower, ci_upper = np.percentile(best_a_mean_dist, 2.5), np.percentile(best_a_mean_dist, 97.5)
        print(f"  - 平均 F1 算术平均数: {mean_a_f1:.4f} (95% CI: [{ci_lower:.4f}, {ci_upper:.4f}])")

        print("\n  - 各基元类型的平均F1分数 (在最佳参数下):")
        f1_scores_by_motif = defaultdict(list)
        for replicate_scores in best_f1_scores_per_replicate:
            for motif_type, f1_score in replicate_scores.items():
                f1_scores_by_motif[motif_type].append(f1_score)

        for motif_type, f1_list in sorted(f1_scores_by_motif.items()):
            mean_f1 = np.mean(f1_list)
            print(f"    - {motif_type.capitalize():<6}: {mean_f1:.4f}")

        if len(results_df_sorted) > 1:
            second_best_row = results_df_sorted.iloc[1]
            second_best_a_mean_dist = second_best_row['a_mean_dist']
            
            print("\n--- 统计显著性分析 (配对检验) ---")
            sec_param_string = f"min-identity={second_best_row['min_identity']}, max-error={int(second_best_row['max_error'])}, min-seed-repeat={int(second_best_row['min_seed_repeat'])}"
            print(f"  - 次优参数组合: {sec_param_string} (平均F1算术平均数: {np.mean(second_best_a_mean_dist):.4f})")

            if np.allclose(best_a_mean_dist, second_best_a_mean_dist):
                print("  - 结论: 无法进行统计比较 (F1分布完全相同)。")
            else:
                res = stats.wilcoxon(best_a_mean_dist, second_best_a_mean_dist, alternative='greater', zero_method='pratt')
                p_value = res.pvalue
                print(f"  - Wilcoxon符号秩检验 p-value (最优 vs. 次优): {p_value:.4f}")
                if p_value < 0.05:
                    print("  - 结论: 最佳参数组合的性能在统计学上显著优于次优组合 (p < 0.05)。")
                else:
                    print("  - 结论: 最佳与次优参数组合的性能无统计学上的显著差异 (p >= 0.05)。")
        else:
            print("\n--- 统计显著性分析 ---", "\n  - 无次优参数可供比较。")

    print("="*80)
    
    motif_f1_cols = [f'f1_{mtype}_avg' for mtype in motifs_to_test.keys()]
    final_cols_order = ['min_identity', 'max_error', 'min_seed_repeat'] + motif_f1_cols + ['f1_mean_avg', 'f1_mean_ci']
    
    cols_to_save = [c for c in final_cols_order if c in results_df.columns]
    results_df[cols_to_save].to_csv(args.output_csv, index=False, float_format="%.4f")
    print(f"\n详细的网格搜索结果已保存到: {args.output_csv}")
    
    if not args.keep_files:
        shutil.rmtree(work_dir)
        print(f"临时工作目录 '{work_dir}' 已被删除。")

if __name__ == "__main__":
    main()
