import unittest

from optimize_graph import compute_counts_by_motif


class GraphSSRBenchmarkScoringTests(unittest.TestCase):
    def test_period_is_required_for_hits_but_all_unmatched_predictions_count(self):
        seq_id = "sim_di_rep0"
        motifs = {"di": "AT"}
        truth = [{"SequenceID": seq_id, "MotifType": "di", "True_Start": 10, "True_End": 29}]
        correct = {"start": 10, "end": 29, "period": 2}
        wrong_period = {"start": 10, "end": 29, "period": 3}
        missing_period = {"start": 10, "end": 29}
        below_iou = {"start": 15, "end": 34, "period": 2}

        cases = (
            ("wrong period only", [wrong_period], truth, (0, 1, 1)),
            ("hit with extra predictions", [wrong_period, correct, missing_period], truth, (1, 2, 0)),
            ("no truth", [wrong_period, missing_period], [], (0, 2, 0)),
            ("low IoU", [below_iou], truth, (0, 1, 1)),
            ("no prediction", [], truth, (0, 0, 1)),
        )
        for label, preds, truths, (tp, fp, fn) in cases:
            with self.subTest(case=label):
                counts = compute_counts_by_motif(
                    {seq_id: preds}, truths, motifs, 0.7, "interval+period"
                )
                self.assertEqual(counts["di"], {"tp": tp, "fp": fp, "fn": fn})

    def test_interval_only_matches_without_a_correct_period(self):
        seq_id = "sim_di_rep0"
        counts = compute_counts_by_motif(
            {seq_id: [{"start": 10, "end": 29, "period": 3}]},
            [{"SequenceID": seq_id, "MotifType": "di", "True_Start": 10, "True_End": 29}],
            {"di": "AT"},
            0.7,
            "interval",
        )
        self.assertEqual(counts["di"], {"tp": 1, "fp": 0, "fn": 0})
