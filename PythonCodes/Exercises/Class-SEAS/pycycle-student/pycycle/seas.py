import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import toms748
from .mesh import num_fault_elements
from .bem import assemble, rhs_op


class FaultMap:
    """Map on-fault indices to bem indices."""
    def __init__(self, mesh):
        """Constructor.

        :param mesh: List of line elements.
        """
        self.Nf = num_fault_elements(mesh)
        self.map = np.ndarray((self.Nf, ), dtype=int)
        f = 0
        for i in range(len(mesh)):
            if mesh[i].is_fault:
                self.map[f] = i
                f = f + 1

    def __len__(self):
        """Number of on-fault elements."""
        return self.Nf

    def __call__(self, f):
        """Get bem index for on-fault index f."""
        return self.map[f]


class IFaultMap:
    """Map bem indices to on-fault indices. 'Inverse' of FaultMap."""
    def __init__(self, mesh):
        """Constructor.

        :param mesh: List of line elements.
        """
        self.N = len(mesh)
        self.imap = np.ndarray((self.N, ), dtype=int)
        f = 0
        for i in range(len(mesh)):
            if mesh[i].is_fault:
                self.imap[i] = f
                f = f + 1
            else:
                self.imap[i] = -1

    def __len__(self):
        """Number of bem elements (= number of line elements in mesh)."""
        return self.N

    def __call__(self, i):
        """Get on-fault index for bem index i."""
        return self.imap[i]


class VariableParams:
    """Prestress and a parameter for rate and state friction.
       May depend on position in space.
    """
    def __init__(self, mesh, a, tau_pre):
        """Constructor.

        :param mesh: List of line elements.
        :param a: Functional R^2 -> R for a parameter.
        :param tau_pre: Functional R^2 -> R for pre-stress.
        """
        Nf = num_fault_elements(mesh)
        self.x = np.ndarray((Nf, 2))
        f = 0
        for i in range(len(mesh)):
            if mesh[i].is_fault:
                self.x[f, :] = mesh[i].collocation_point()
                f = f + 1
        self.a = np.ndarray((Nf, ))
        self.tau_pre = np.ndarray((Nf, ))
        for f in range(Nf):
            self.a[f] = a(self.x[f, :])
            self.tau_pre[f] = tau_pre(self.x[f, :])


class ConstantParams:
    """Constant parameters for rate and state friction."""
    def __init__(self, rho, v_s, Vp, V0, b, L, f0, sn, Vinit):
        """Constructor.

        :param rho: density [g/m^3]
        :param v_s: shear wave velocity [km/s]
        :param Vp: pate rate [m/s]
        :param V0: reference slip rate [m/s]
        :param b: b parameter
        :param L: critical slip distance (aka D_c) [m]
        :param f0: reference friction coefficient
        :param sn: normal stress [MPa]
        :param Vinit: initial slip rate [m/s]
        """
        self.mu = rho * v_s**2
        self.eta = rho * v_s / 2.0
        self.Vp = Vp
        self.V0 = V0
        self.b = b
        self.L = L
        self.f0 = f0
        self.sn = sn
        self.Vinit = Vinit


class Context:
    """Context which contains everything we need in evaluating the right-hand side
       in the time integrator.
    """
    def __init__(self, mesh, G, dG_dn, vp, cp):
        """Constructor.

        :param mesh: List of line elements 
        :param G: Green's function
        :param dG_dn: Directional derivative of Green's function
        :param vp: VariableParams
        :param cp: ConstantParams
        """
        self.A = assemble(G, mesh)
        self.B = rhs_op(dG_dn, mesh)
        self.lu, self.piv = lu_factor(self.A)
        self.map = FaultMap(mesh)
        self.imap = IFaultMap(mesh)
        self.vp = vp
        self.cp = cp

    def traction(self, time, u):
        """Computes tau at time 'time' for on-fault displacement u (u has size Nf).

        :param time: Time [s]
        :param u: Displacement vector [m]
        """
        N = len(self.imap)
        Nf = len(self.map)
        tau = np.zeros((Nf, ))
        g = np.zeros((N, ))

        # TODO: implement

        return self.vp.tau_pre + self.cp.mu * tau

    def psi0(self, f):
        """Compute initial state for the f-th on-fault element.

        :param f: Fault element number
        """
        tau = -self.vp.tau_pre[f]
        a = self.vp.a[f]
        s = np.sinh((tau - self.cp.eta * self.cp.Vinit) / (a * self.cp.sn))
        l = np.log((2.0 * self.cp.V0 / self.cp.Vinit) * s)
        return a * l

    def friction_law(self, f, V, psi):
        """Evaluate friction law.

        :param f: Fault element number
        :param V: Slip-rate (scalar)
        :param psi: State (scalar)
        """
        a = self.vp.a[f]
        e = np.exp(psi / a)
        f = a * np.arcsinh((V / (2.0 * self.cp.V0)) * e)
        return self.cp.sn * f

    def slip_rate(self, f, tau, psi):
        """Obtain slip-rate by solving tau + friction_law(V, psi) + eta V = 0 for V.

        :param f: Fault element number
        :param tau: Traction (scalar)
        :param psi: State (scalar)
        """
        # TODO: implement
        return 0.0

    def state_law(self, f, V, psi):
        """Evaluate ageing law.

           psi  = f0 + b log(V0 theta / L)
             -->  L / V0 exp( ( psi - f0 ) / b ) = theta
           psi' = b theta' / theta = b / theta - b V / L
                = b V0 / L exp( ( f0 - psi ) / b ) - b V / L
                = b V0 / L ( exp( ( f0 - psi ) / b ) - V / V0 )

        :param f: Fault element number
        :param V: Slip-rate (scalar)
        :param psi: State (scalar)
        """
        return self.cp.b * self.cp.V0 / self.cp.L * (np.exp(
            (self.cp.f0 - psi) / self.cp.b) - V / self.cp.V0)


def y0(ctx):
    """Compute initial condition for parameters set in context.
       Note: Initial displacement is zero.

    :param ctx: Context
    """
    Nf = len(ctx.map)
    y0 = np.zeros((2 * Nf, ))
    for f in range(Nf):
        y0[2 * f + 1] = ctx.psi0(f)
    return y0


def F(t, y, ctx, callback=None):
    """Evaluate right-hand side of the SEAS ODE. The state vector y interleaves displacement
       and psi variable in the following way:
       [S_0, psi_0, S_1, psi_1, ..., S_{N_f}, psi_{N_f}],
       where N_f is the number of fault elements.

    :param t: Time
    :param y: State vector according to solve_ivp (not to be confused with state variable psi).
    :param ctx: Context
    :param callback: Callback with signature (t, S, V, psi, tau)
    """
    tau = ctx.traction(t, y[::2])
    Nf = y.shape[0] // 2
    fy = np.ndarray(y.shape)
    for f in range(Nf):
        psi = y[2 * f + 1]
        V = ctx.slip_rate(f, tau[f], psi)
        fy[2 * f] = V
        fy[2 * f + 1] = ctx.state_law(f, V, psi)
    if callback is not None:
        callback(t,y[::2],fy[::2],y[1::2],tau)
    return fy
