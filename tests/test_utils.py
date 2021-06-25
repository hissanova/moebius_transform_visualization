from unittest import TestCase, main

from numpy import inf, pi

from src.utils import cot


class CotTestCases(TestCase):
    def test_positive_infinity(self):
        res = cot(0)
        self.assertEqual(res, +inf)
        self.assertTrue(res > 0)

    def test_negative_infinity(self):
        res = cot(pi)
        self.assertEqual(res, -inf)
        self.assertTrue(res < 0)

    def test_exceptions(self):
        self.assertRaises(Exception, cot, -1)
        self.assertRaises(Exception, cot, 2 * pi)


if __name__ == '__main__':
    main()
