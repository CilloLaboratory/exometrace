from pathlib import Path
import tempfile
import unittest

from scripts.validate_samplesheet import read_samplesheet, validate_rows


class SampleSheetTests(unittest.TestCase):
    def test_valid_samplesheet_passes(self) -> None:
        path = Path("tests/test_data/samplesheet_valid.csv")
        rows, _ = read_samplesheet(path)
        self.assertEqual(validate_rows(rows, path), [])

    def test_invalid_samplesheet_fails(self) -> None:
        path = Path("tests/test_data/samplesheet_invalid.csv")
        rows, _ = read_samplesheet(path)
        errors = validate_rows(rows, path)
        self.assertGreaterEqual(len(errors), 3)

    def test_missing_file_detection(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            sheet = Path(tempdir) / "samples.csv"
            sheet.write_text(
                "patient_id,tumor_sample,normal_sample,tumor_r1,tumor_r2,normal_r1,normal_r2,bait_bed\n"
                "P001,T,N,a.fastq.gz,b.fastq.gz,c.fastq.gz,d.fastq.gz,targets.bed\n",
                encoding="utf-8",
            )
            rows, _ = read_samplesheet(sheet)
            errors = validate_rows(rows, sheet, check_files=True)
            self.assertTrue(any("missing FASTQ" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
