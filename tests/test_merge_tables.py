import csv
import tempfile
import unittest

from pathlib import Path

from scripts.merge_tables import merge_identical_headers, merge_union_by_first_column


class MergeTablesTests(unittest.TestCase):
    def test_merge_identical_headers_rejects_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            first = root / "a.tsv"
            second = root / "b.tsv"
            output = root / "out.tsv"
            first.write_text("patient_id\tTP53\nP1\t1\n", encoding="utf-8")
            second.write_text("patient_id\tKRAS\nP2\t1\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                merge_identical_headers([first, second], output)

    def test_merge_union_by_first_column_fills_missing_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            first = root / "a.tsv"
            second = root / "b.tsv"
            output = root / "out.tsv"
            first.write_text("patient_id\tTP53\nT2D_4605_POST\t1\n", encoding="utf-8")
            second.write_text("patient_id\tKRAS\nT2D_4605_PRE\t1\n", encoding="utf-8")

            merge_union_by_first_column([first, second], output, "0")

            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["patient_id"], "T2D_4605_POST")
            self.assertEqual(rows[0]["TP53"], "1")
            self.assertEqual(rows[0]["KRAS"], "0")
            self.assertEqual(rows[1]["patient_id"], "T2D_4605_PRE")
            self.assertEqual(rows[1]["TP53"], "0")
            self.assertEqual(rows[1]["KRAS"], "1")


if __name__ == "__main__":
    unittest.main()
