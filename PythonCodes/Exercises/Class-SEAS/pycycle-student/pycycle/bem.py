import numpy as np
from scipy.integrate import quad


def assemble(G, mesh):
    """Assemble the BEM operator A (left-hand side).

    :param G: Green's function G(x, xi)
    :param mesh: List of line elements
    """
    M = len(mesh)
    A = np.ndarray((M, M))

    # TODO: implement

    return A


def rhs_op(dG_dn, mesh):
    """Assemble the BEM operator B, which is used to construct the right-hand side,
       i.e. b = B @ u.

    :param dG_dn: Directional derivative of Green's function dG_dn(x, xi, n)
    :param mesh: List of line elements
    """
    M = len(mesh)
    B = 0.5 * np.eye(M)
    for i in range(M):
        xc = mesh[i].collocation_point()
        for j in range(M):
            K = lambda t: dG_dn(xc, mesh[j].xi(t), mesh[j].n) * mesh[j].factor(
                t)
            B[i, j] += quad(K, -1, 0)[0] + quad(K, 0, 1)[0]
    return B
