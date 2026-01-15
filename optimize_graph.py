#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import argparse
import random
import textwrap
import shutil
import json
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    print("Error: numpy is required. Please run: pip install numpy")
    raise SystemExit(1)

try:
    from scipy import stats
except ImportError:
    stats = None


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


def load_truth_data(truth_file):
    records = []
    try:
        with open(truth_file, 'r') as f:
            header = f.readline().strip().split('\t')
            if not header or "SequenceID" not in header:
                return []
            for line in f:
                parts = line.rstrip("\n").split('\t')
                if len(parts) < 4:
                    continue
                try:
                    records.append({
                        'SequenceID': parts[0],
                        'MotifType': parts[1],
                        'True_Start': int(parts[2]),
                        'True_End': int(parts[3]),
                    })
                except ValueError:
                    continue
    except OSError:
        return []
    return records


def parse_ssrs_tsv(ssrs_path):
    predictions_by_seq = defaultdict(list)
    if not os.path.exists(ssrs_path) or os.path.getsize(ssrs_path) == 0:
        return predictions_by_seq
    try:
        with open(ssrs_path, 'r') as f:
            header = f.readline().strip().split('\t')
            if not header:
                return predictions_by_seq
            try:
                idx_seq = header.index("Sequence_ID")
                idx_start = header.index("Start")
                idx_end = header.index("End")
                idx_motif = header.index("Inferred_Motif") if "Inferred_Motif" in header else None
            except ValueError:
                return predictions_by_seq
            for line in f:
                parts = line.rstrip("\n").split('\t')
                if len(parts) <= max(idx_seq, idx_start, idx_end):
                    continue
                seq_id = parts[idx_seq]
                try:
                    start = int(parts[idx_start])
                    end = int(parts[idx_end])
                except ValueError:
                    continue
                period = None
                if idx_motif is not None and len(parts) > idx_motif:
                    motif = parts[idx_motif].strip()
                    if motif:
                        period = len(motif)
                predictions_by_seq[seq_id].append({'start': start, 'end': end, 'period': period})
    except OSError:
        return defaultdict(list)
    return predictions_by_seq


def compute_counts_by_motif(predictions_by_seq, truth_data, motifs_to_test, overlap_threshold, match_mode):
    counts_by_motif = {mtype: {'tp': 0, 'fp': 0, 'fn': 0} for mtype in motifs_to_test}
    truths_by_motif = defaultdict(list)
    for t in truth_data:
        truths_by_motif[t['MotifType']].append(t)
    motif_periods = {mtype: len(motif) for mtype, motif in motifs_to_test.items()}

    for motif_type in motifs_to_test:
        tp = fp = fn = 0
        truths_for_motif = truths_by_motif.get(motif_type, [])
        truth_by_seq = {t['SequenceID']: t for t in truths_for_motif}
        truth_seq_ids = set(truth_by_seq.keys())
        pred_seq_ids_for_motif = {sid for sid in predictions_by_seq if motif_type in sid}
        all_seq_ids = truth_seq_ids.union(pred_seq_ids_for_motif)
        required_period = motif_periods[motif_type] if match_mode == "interval+period" else None

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
                    fp += max(0, sum(1 for p in preds if required_period is None or p.get('period') == required_period) - 1)
                else:
                    fn += 1
                    fp += sum(1 for p in preds if required_period is None or p.get('period') == required_period)
            else:
                fp += sum(1 for p in preds if required_period is None or p.get('period') == required_period)

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


def compute_null_overlaps(predictions_by_seq, truth_by_seq, seq_lengths, rng, motifs_to_test, match_mode):
    overlaps = []
    motif_periods = {mtype: len(motif) for mtype, motif in motifs_to_test.items()}
    for seq_id, preds in predictions_by_seq.items():
        truth = truth_by_seq.get(seq_id)
        seq_len = seq_lengths.get(seq_id)
        if not truth or not seq_len:
            continue
        required_period = None
        if match_mode == "interval+period":
            required_period = motif_periods.get(truth.get('MotifType'))
        true_start, true_end = truth['True_Start'], truth['True_End']
        for pred in preds:
            if required_period is not None and pred.get('period') != required_period:
                continue
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


def run_command(cmd, cwd=None, env=None):
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd, text=True, env=env)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        raise RuntimeError(f"Command failed: {' '.join(cmd)}; stderr={stderr}")


def build_graph_command(graph_path, fasta_path, output_prefix, params):
    base_cmd = [
        graph_path, fasta_path, output_prefix,
        "-k", str(params['kmer_size']),
        "--min_len", str(params['min_len']),
        "--min_cov", str(params['min_cov']),
        "--max_merge_dist", str(params['max_merge_dist']),
    ]
    if graph_path.endswith(".py"):
        return ["python3"] + base_cmd
    if os.access(graph_path, os.X_OK):
        return base_cmd
    raise RuntimeError(f"graph_path is not executable: {graph_path}")


def run_ssr_graph(graph_path, fasta_path, output_prefix, params, env):
    cmd = build_graph_command(graph_path, fasta_path, output_prefix, params)
    run_command(cmd, env=env)
    ssrs_path = f"{output_prefix}.ssrs"
    return parse_ssrs_tsv(ssrs_path)


def parse_range_str(input_str):
    return [p.strip() for p in input_str.split(',') if p.strip()]


def parse_range_int(input_str):
    return [int(p) for p in parse_range_str(input_str)]


def parse_range_float(input_str):
    return [float(p) for p in parse_range_str(input_str)]


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


def main():
    parser = argparse.ArgumentParser(
        description="Grid search and evaluation for ssr_graph.py",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--graph_path", default="/app/wtSSR/ssr_graph.py")
    parser.add_argument("--num_replicates", type=int, default=30)
    parser.add_argument("--snp_rate", type=float, default=0.02)
    parser.add_argument("--indel_rate", type=float, default=0.01)
    parser.add_argument("--flank_len", type=int, default=200)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output_csv", default="ssr_graph_grid_search_results.csv")
    parser.add_argument("--keep_files", action='store_true')
    parser.add_argument("--reuse_sim", action='store_true')
    parser.add_argument("--force_regen", action='store_true')
    parser.add_argument("--max_workers", type=int, default=1)
    parser.add_argument("--overlap_threshold", type=float, default=0.5)
    parser.add_argument("--match_mode", choices=["interval", "interval+period"], default="interval",
                        help="评估模式: interval=仅区间重叠; interval+period=区间重叠且period匹配。")
    parser.add_argument("--calibrate_threshold", action='store_true')
    parser.add_argument("--calibration_fraction", type=float, default=0.2)
    parser.add_argument("--calibration_replicates", type=int, default=None)
    parser.add_argument("--threshold_grid", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    parser.add_argument("--target_precision", type=float, default=0.9)
    parser.add_argument("--null_quantile", type=float, default=0.95)
    parser.add_argument("--disable_pyfastx", action='store_true')

    parser.add_argument("--graph_k_range", default="21")
    parser.add_argument("--graph_min_len_range", default="30")
    parser.add_argument("--graph_min_cov_range", default="2.0")
    parser.add_argument("--graph_max_merge_dist_range", default="25")
    args = parser.parse_args()

    if args.reuse_sim and args.force_regen:
        print("Error: --reuse_sim and --force_regen cannot be used together.")
        return
    if args.max_workers < 1:
        print("Error: --max_workers must be >= 1.")
        return
    if not (0 < args.overlap_threshold <= 1):
        print("Error: --overlap_threshold must be in (0, 1].")
        return
    if not os.path.exists(args.graph_path):
        print(f"Error: script not found: {args.graph_path}")
        return
    if not args.graph_path.endswith(".py") and not os.access(args.graph_path, os.X_OK):
        print(f"Error: graph_path is not executable: {args.graph_path}")
        return

    threshold_grid = parse_threshold_grid(args.threshold_grid)
    if not threshold_grid or any(t <= 0 or t > 1 for t in threshold_grid):
        print("Error: --threshold_grid must be in (0, 1].")
        return

    if args.calibrate_threshold:
        if args.calibration_replicates is not None:
            calibration_count = args.calibration_replicates
        else:
            calibration_count = max(1, int(args.num_replicates * args.calibration_fraction))
        if calibration_count >= args.num_replicates:
            print("Error: calibration set must be smaller than total replicates.")
            return
        calibration_indices = list(range(calibration_count))
        evaluation_indices = list(range(calibration_count, args.num_replicates))
        print(f"Calibration replicates: {len(calibration_indices)}; evaluation replicates: {len(evaluation_indices)}")
    else:
        calibration_indices = []
        evaluation_indices = list(range(args.num_replicates))

    motifs_to_test = {'mono': 'A', 'di': 'AT', 'tri': 'AGC', 'tetra': 'AATC', 'penta': 'AGCTC', 'hexa': 'AATGCC'}
    repeats_map = {'mono': 50, 'di': 25, 'tri': 17, 'tetra': 12, 'penta': 10, 'hexa': 8}
    work_dir = "ssr_graph_workspace"
    os.makedirs(work_dir, exist_ok=True)

    env = os.environ.copy()
    if args.disable_pyfastx or args.max_workers > 1:
        env["PYFASTX_DISABLE"] = "1"

    def prepare_sim_data(replicate_id):
        fasta_file = os.path.join(work_dir, f"combined_rep_{replicate_id}.fasta")
        truth_file = os.path.join(work_dir, f"combined_truth_rep_{replicate_id}.tsv")
        meta_file = os.path.join(work_dir, f"combined_rep_{replicate_id}.meta.json")

        replicate_seed = args.seed + replicate_id if args.seed is not None else None
        sim_meta = build_sim_metadata(
            motifs_to_test, repeats_map, args.flank_len,
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

        if not reuse_sim:
            if replicate_seed is not None:
                random.seed(replicate_seed)
                np.random.seed(replicate_seed)
            generate_combined_dataset(
                motifs_to_test, repeats_map, args.flank_len,
                args.snp_rate, args.indel_rate, work_dir, replicate_id=replicate_id
            )
            write_sim_metadata(meta_file, sim_meta)

        seq_lengths = read_fasta_lengths(fasta_file)
        return fasta_file, truth_file, seq_lengths

    graph_combos = []
    for k in parse_range_int(args.graph_k_range):
        for min_len in parse_range_int(args.graph_min_len_range):
            for min_cov in parse_range_float(args.graph_min_cov_range):
                for max_merge_dist in parse_range_int(args.graph_max_merge_dist_range):
                    graph_combos.append({
                        'kmer_size': k,
                        'min_len': min_len,
                        'min_cov': min_cov,
                        'max_merge_dist': max_merge_dist
                    })

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

        print("\n--- Calibration ---")
        for idx, i in enumerate(calibration_indices, start=1):
            print(f"\n--- Calibration Replicate {idx}/{len(calibration_indices)} (original {i + 1}) ---")
            fasta_file, truth_file, seq_lengths = prepare_sim_data(i)
            truth_data = load_truth_data(truth_file)
            truth_by_seq = {t['SequenceID']: t for t in truth_data}

            def run_combo_calibration(j, params):
                output_prefix = os.path.join(work_dir, f"graph_combo_{j}_rep_{i}")
                predictions_by_seq = run_ssr_graph(args.graph_path, fasta_file, output_prefix, params, env)
                counts_by_threshold = {}
                for t in threshold_grid:
                    counts_by_threshold[t] = compute_counts_by_motif(
                        predictions_by_seq,
                        truth_data,
                        motifs_to_test,
                        t,
                        args.match_mode
                    )
                rng = random.Random(null_seed_base + i * 100000 + j)
                null_ov = compute_null_overlaps(
                    predictions_by_seq,
                    truth_by_seq,
                    seq_lengths,
                    rng,
                    motifs_to_test,
                    args.match_mode
                )
                return counts_by_threshold, null_ov

            if args.max_workers == 1:
                for j, params in enumerate(graph_combos):
                    print(f"\r    - Testing combo {j + 1}/{len(graph_combos)}...", end="")
                    counts_by_threshold, null_ov = run_combo_calibration(j, params)
                    for t, counts in counts_by_threshold.items():
                        for motif_type, c in counts.items():
                            total = threshold_stats[t][motif_type]
                            total['tp'] += c['tp']
                            total['fp'] += c['fp']
                            total['fn'] += c['fn']
                    null_overlaps.extend(null_ov)
                print("\r    - Completed all combos.           ")
            else:
                completed = 0
                with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                    futures = {
                        executor.submit(run_combo_calibration, j, params): j
                        for j, params in enumerate(graph_combos)
                    }
                    for future in as_completed(futures):
                        counts_by_threshold, null_ov = future.result()
                        for t, counts in counts_by_threshold.items():
                            for motif_type, c in counts.items():
                                total = threshold_stats[t][motif_type]
                                total['tp'] += c['tp']
                                total['fp'] += c['fp']
                                total['fn'] += c['fn']
                        null_overlaps.extend(null_ov)
                        completed += 1
                        print(f"\r    - Testing combo {completed}/{len(graph_combos)}...", end="")
                print("\r    - Completed all combos.           ")

        if null_overlaps:
            t_null = float(np.quantile(null_overlaps, args.null_quantile))
        else:
            t_null = 0.0

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
        with open(curve_path, 'w') as f:
            f.write("threshold,macro_precision,macro_recall,macro_f1\n")
            for row in threshold_records:
                f.write(f"{row['threshold']:.4f},{row['macro_precision']:.4f},{row['macro_recall']:.4f},{row['macro_f1']:.4f}\n")
        meta = {
            't_null': t_null,
            'selected_threshold': selected_threshold,
            'threshold_grid': threshold_grid,
            'target_precision': args.target_precision,
            'null_quantile': args.null_quantile,
            'calibration_replicates': len(calibration_indices),
            'evaluation_replicates': len(evaluation_indices),
            'match_mode': args.match_mode,
        }
        with open(meta_path, 'w') as f:
            json.dump(meta, f, ensure_ascii=True, indent=2, sort_keys=True)

    print("\n--- Evaluation ---")
    for idx, i in enumerate(evaluation_indices, start=1):
        print(f"\n--- Evaluation Replicate {idx}/{len(evaluation_indices)} (original {i + 1}) ---")
        fasta_file, truth_file, _ = prepare_sim_data(i)
        truth_data = load_truth_data(truth_file)

        def run_combo_eval(j, params):
            output_prefix = os.path.join(work_dir, f"graph_combo_{j}_rep_{i}")
            predictions_by_seq = run_ssr_graph(args.graph_path, fasta_file, output_prefix, params, env)
            counts_by_motif = compute_counts_by_motif(
                predictions_by_seq,
                truth_data,
                motifs_to_test,
                selected_threshold,
                args.match_mode
            )
            f1_scores, macro_precision, macro_recall, macro_f1 = compute_metrics_from_counts(
                counts_by_motif, motifs_to_test
            )
            return j, f1_scores, macro_precision, macro_recall, macro_f1

        if args.max_workers == 1:
            for j, params in enumerate(graph_combos):
                print(f"\r    - Testing combo {j + 1}/{len(graph_combos)}...", end="")
                j_idx, f1_scores, macro_precision, macro_recall, macro_f1 = run_combo_eval(j, params)
                results_dist[j_idx].append({
                    'f1_scores': f1_scores,
                    'macro_precision': macro_precision,
                    'macro_recall': macro_recall,
                    'macro_f1': macro_f1,
                })
            print("\r    - Completed all combos.           ")
        else:
            completed = 0
            with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                futures = {
                    executor.submit(run_combo_eval, j, params): j
                    for j, params in enumerate(graph_combos)
                }
                for future in as_completed(futures):
                    j_idx, f1_scores, macro_precision, macro_recall, macro_f1 = future.result()
                    results_dist[j_idx].append({
                        'f1_scores': f1_scores,
                        'macro_precision': macro_precision,
                        'macro_recall': macro_recall,
                        'macro_f1': macro_f1,
                    })
                    completed += 1
                    print(f"\r    - Testing combo {completed}/{len(graph_combos)}...", end="")
            print("\r    - Completed all combos.           ")

    all_rows = []
    summary_rows = []
    a_mean_dist_map = {}
    for j, params in enumerate(graph_combos):
        records = results_dist.get(j, [])
        macro_f1_list = [r['macro_f1'] for r in records]
        a_mean_dist_map[j] = macro_f1_list
        if macro_f1_list:
            avg_macro_f1 = float(np.mean(macro_f1_list))
            ci_lower = float(np.percentile(macro_f1_list, 2.5))
            ci_upper = float(np.percentile(macro_f1_list, 97.5))
            ci_str = f"[{ci_lower:.4f}, {ci_upper:.4f}]"
        else:
            avg_macro_f1 = 0.0
            ci_str = "[0.0000, 0.0000]"

        f1_scores_by_motif = defaultdict(list)
        for rec in records:
            for motif_type, f1_score in rec['f1_scores'].items():
                f1_scores_by_motif[motif_type].append(f1_score)

        row = {
            'kmer_size': params['kmer_size'],
            'min_len': params['min_len'],
            'min_cov': params['min_cov'],
            'max_merge_dist': params['max_merge_dist'],
            'f1_mean_avg': avg_macro_f1,
            'f1_mean_ci': ci_str,
        }
        for motif_type in motifs_to_test:
            f1_list = f1_scores_by_motif.get(motif_type, [])
            row[f'f1_{motif_type}_avg'] = float(np.mean(f1_list)) if f1_list else 0.0
        all_rows.append(row)
        summary_rows.append((avg_macro_f1, j, row))

    all_rows_sorted = sorted(all_rows, key=lambda r: r.get('f1_mean_avg', 0), reverse=True)
    if not all_rows_sorted:
        print("No results to report.")
        return

    print("\n--- Top Results ---")
    for row in all_rows_sorted[:5]:
        print(row)

    summary_rows = sorted(summary_rows, key=lambda item: item[0], reverse=True)
    best_avg, best_idx, best_row = summary_rows[0]
    second = summary_rows[1] if len(summary_rows) > 1 else None

    print("\n" + "=" * 80)
    print("--- 最终结论: 推荐使用的最佳统一Graph参数组合 ---")
    print("  - 最佳参数组合:")
    print(f"    - -k: {best_row['kmer_size']}")
    print(f"    - --min_len: {best_row['min_len']}")
    print(f"    - --min_cov: {best_row['min_cov']}")
    print(f"    - --max_merge_dist: {best_row['max_merge_dist']}")
    print(f"  - 平均 F1 算术平均数: {best_row['f1_mean_avg']:.4f} (95% CI: {best_row['f1_mean_ci']})")
    print("\n  - 各基元类型的平均F1分数 (在最佳参数下):")
    for motif_type in motifs_to_test:
        label = motif_type.capitalize()
        print(f"    - {label:<6}: {best_row.get(f'f1_{motif_type}_avg', 0.0):.4f}")

    if second is not None:
        second_avg, second_idx, second_row = second
        print("\n--- 统计显著性分析 (配对检验) ---")
        second_param_str = (
            f"-k={second_row['kmer_size']}, "
            f"--min_len={second_row['min_len']}, "
            f"--min_cov={second_row['min_cov']}, "
            f"--max_merge_dist={second_row['max_merge_dist']}"
        )
        print(f"  - 次优参数组合: {second_param_str} (平均F1算术平均数: {second_avg:.4f})")

        best_dist = a_mean_dist_map.get(best_idx, [])
        second_dist = a_mean_dist_map.get(second_idx, [])
        if not best_dist or not second_dist:
            print("  - 结论: 无法进行统计比较 (F1分布为空)。")
        elif stats is None:
            print("  - 结论: 未安装SciPy，跳过统计检验。")
        elif np.allclose(best_dist, second_dist):
            print("  - 结论: 无法进行统计比较 (F1分布完全相同)。")
        else:
            try:
                res = stats.wilcoxon(best_dist, second_dist, alternative='greater', zero_method='pratt')
                p_value = res.pvalue
                print(f"  - Wilcoxon符号秩检验 p-value (最优 vs. 次优): {p_value:.4f}")
                if p_value < 0.05:
                    print("  - 结论: 最佳参数组合的性能在统计学上显著优于次优组合 (p < 0.05)。")
                else:
                    print("  - 结论: 最佳与次优参数组合的性能无统计学上的显著差异 (p >= 0.05)。")
            except ValueError:
                print("  - 结论: 无法进行统计比较 (配对差异全为0)。")
    else:
        print("\n--- 统计显著性分析 ---")
        print("  - 无次优参数可供比较。")
    print("=" * 80)

    fieldnames = [
        'kmer_size', 'min_len', 'min_cov', 'max_merge_dist',
        'f1_mono_avg', 'f1_di_avg', 'f1_tri_avg', 'f1_tetra_avg', 'f1_penta_avg', 'f1_hexa_avg',
        'f1_mean_avg', 'f1_mean_ci'
    ]
    with open(args.output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(fieldnames)
        for row in all_rows_sorted:
            values = [row.get(k, "") for k in fieldnames]
            writer.writerow(values)

    print(f"\nSaved results to: {args.output_csv}")
    if not args.keep_files:
        shutil.rmtree(work_dir)
        print(f"Removed workspace: {work_dir}")


if __name__ == "__main__":
    main()
