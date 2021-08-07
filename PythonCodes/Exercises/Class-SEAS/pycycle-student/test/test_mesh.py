import numpy as np
import unittest
from scipy.integrate import quad

from pycycle.mesh import LineElement, InfiniteLineElement


class TestLineElement(unittest.TestCase):
    def setUp(self):
        self.a = np.array((1, 1))
        self.b = np.array((3, 3))
        self.midpoint = np.array((2, 2))
        self.normal = (-1, 1)
        self.line = LineElement(self.a, self.b, self.normal, False)

    def test_xi(self):
        self.assertTrue(all(self.a == self.line.xi(-1)))
        self.assertTrue(all(self.b == self.line.xi(1)))
        self.assertTrue(all(self.midpoint == self.line.xi(0)))

    def test_basis(self):
        self.assertEqual(self.line.basis(-1), 1)
        self.assertEqual(self.line.basis(0), 1)
        self.assertEqual(self.line.basis(1), 1)

    def test_integral(self):
        """Tests int_{\Gamma} ln|x| dx, where Gamma is the
           line from self.a to self.b.
        """

        K = lambda t: np.log(np.linalg.norm(self.line.xi(t))
                             ) * self.line.factor(t)
        a_norm = np.linalg.norm(self.a)
        analytic = (3 * (np.log(3) - 1) - (-1) +
                    (3 - 1) * np.log(a_norm)) * a_norm
        numeric = quad(K, -1, 1)[0]
        self.assertAlmostEqual(analytic, numeric)


class TestInfiniteLineElement(unittest.TestCase):
    def setUp(self):
        self.a = (1, 1)
        self.normal = (-1, 1)
        self.line = InfiniteLineElement(self.a, self.normal)

    def test_xi(self):
        self.assertTrue(all(self.a == self.line.xi(-1)))
        xi0 = self.line.xi(0)
        for i in range(2):
            self.assertAlmostEqual(3., xi0[i])
        xi0_99 = self.line.xi(0.99)
        for i in range(2):
            self.assertAlmostEqual(399., xi0_99[i])

    def test_basis(self):
        self.assertAlmostEqual(self.line.basis(-1), 1)
        self.assertAlmostEqual(self.line.basis(0), 1. / 3**2)
        self.assertAlmostEqual(self.line.basis(0.99), 1. / 399.**2)
        self.assertAlmostEqual(self.line.basis(1), 0.0)

    def test_integral(self):
        """Tests int_{\Gamma} ln|x| * |a|**2/|x|**2 dx, where Gamma is
           the line from self.a to infinity (in the direction of self.a).
        """

        K = lambda t: np.log(np.linalg.norm(self.line.xi(t))
                             ) * self.line.factor(t)
        a_norm = np.linalg.norm(self.a)
        analytic = (1 + np.log(a_norm)) * a_norm
        numeric = quad(K, -1, 1)[0]
        self.assertAlmostEqual(analytic, numeric)


if __name__ == '__main__':
    unittest.main()
