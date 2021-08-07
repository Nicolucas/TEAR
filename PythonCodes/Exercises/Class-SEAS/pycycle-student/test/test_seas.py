import numpy as np
import random
import unittest

import pycycle.green as green
from pycycle.mesh import num_fault_elements, LineElement, InfiniteLineElement, tessellate_line
from pycycle.seas import Context, ConstantParams, VariableParams


class TestBEM(unittest.TestCase):
    def setUp(self):
        a = np.array((0, 0.1))
        b = np.array((0, 1))
        normal = (-1, 0)
        self.mesh = [LineElement((0, 0), a, normal, False)]
        self.mesh += tessellate_line(a, b, 0.1, normal, True)
        self.mesh += [InfiniteLineElement(b, normal)]

        rho = 2.670  # density [g/m^3]
        v_s = 3.464  # shear wave velocity [km/s]
        Vp = 1e-9  # plate rate [m/s]
        V0 = 1e-6  # reference slip rate [m/s]
        b = 0.015  # b parameter
        L = 0.014  # critical slip distance [m]
        f0 = 0.6  # reference friction coefficient
        sn = 50  # normal stress [MPa]
        Vinit = 1e-9  # initial slip rate [m/s]
        cp = ConstantParams(rho, v_s, Vp, V0, b, L, f0, sn, Vinit)

        a = lambda x: 0.10
        tau_pre = lambda x: -20
        vp = VariableParams(self.mesh, a, tau_pre)

        self.ctx = Context(self.mesh, green.G_fs, green.dG_fs_dn, vp, cp)

    def test_slip_rate(self):
        def test_constraint(tau, psi):
            V = self.ctx.slip_rate(0, tau, psi)
            return tau + self.ctx.friction_law(0, V, psi) + self.ctx.cp.eta * V

        # Corner cases
        self.assertAlmostEqual(test_constraint(0, 0), 0.0)
        self.assertAlmostEqual(test_constraint(0, 1), 0.0)
        self.assertAlmostEqual(test_constraint(-100, 0), 0.0)
        self.assertAlmostEqual(test_constraint(-100, 1), 0.0)

        # Random trials
        ntrials = 100
        for i in range(ntrials):
            tau = random.uniform(-40, 0)
            psi = random.uniform(0, 1)
            self.assertAlmostEqual(test_constraint(tau, psi), 0.0)

    def test_traction(self):
        M = num_fault_elements(self.mesh)
        u = np.zeros((M, ))
        tau0 = self.ctx.traction(0, u)
        self.assertEqual(tau0.shape[0], M)
        for i in range(M):
            self.assertAlmostEqual(tau0[i], -20.0)

        tau1 = self.ctx.traction(1e6, u)
        tau1_ref = [
            -20.2770392, -20.18965558, -20.18867066, -20.19184425,
            -20.20018175, -20.21425668, -20.23834119, -20.2767227, -20.51650913
        ]
        for i in range(M):
            self.assertAlmostEqual(tau1[i], tau1_ref[i])

        for i in range(M):
            c = (M - 1) / 2
            x = i - c
            u[i] = -np.exp(-c**2 / (c**2 - x**2)) if np.abs(x) < c else 0
        tau2 = self.ctx.traction(0, u)
        tau2_ref = [
            -2.89302998, -13.06762148, -37.0195062, -40.08661435, -40.73051422,
            -40.50087324, -37.88783157, -14.51133166, -4.3230266
        ]
        for i in range(M):
            self.assertAlmostEqual(tau2[i], tau2_ref[i])


if __name__ == '__main__':
    unittest.main()
