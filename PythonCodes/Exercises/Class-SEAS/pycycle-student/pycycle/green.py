import numpy as np


def G(x, xi):
    """Fundamental solution (homogeneous whole-space). """
    return -np.log(np.linalg.norm(x - xi)) / (2 * np.pi)


def dG_dn(x, xi, n):
    """Directional derivative of fundamental solution. """
    d = x - xi
    eps = np.finfo(d.dtype).eps
    if np.abs(np.inner(d, n)) < eps:
        return 0
    else:
        return np.inner(d, n) / (np.inner(d, d) * 2 * np.pi)


def G_fs(x, xi):
    """Green's function for half-space with free surface."""
    x_tilde = x.copy()
    x_tilde[1] = -x_tilde[1]
    return G(x, xi) + G(x_tilde, xi)


def dG_fs_dn(x, xi, n):
    """Directional derivative of G_fs."""
    x_tilde = x.copy()
    x_tilde[1] = -x_tilde[1]
    return dG_dn(x, xi, n) + dG_dn(x_tilde, xi, n)
