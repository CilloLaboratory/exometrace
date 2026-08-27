import unittest

from scripts.calculate_tmb import calculate_tmb
from scripts.calculate_vaf import calculate_vaf


class MetricTests(unittest.TestCase):
    def test_calculate_vaf(self) -> None:
        self.assertAlmostEqual(calculate_vaf(30, 10), 0.25)

    def test_calculate_tmb(self) -> None:
        self.assertAlmostEqual(calculate_tmb(120, 30.0), 4.0)


if __name__ == "__main__":
    unittest.main()
