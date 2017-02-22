"""
Contains the Legendre Polynomial calculation for Schrodinger's Equation.
The return values of each function will return a symbolic equation which can be
lambdified with Sympy to generate a solution. The differential equations will
also be removed for posterity.
"""

from functions import fact
from scipy import special
import numpy as np
import sympy as sp

def c_assoc_legendre(m, l):
    """returns cheat asscociated legendre."""
    return lambda theta: abs(special.lpmn(m, l, np.cos(theta))[0][m][l])
    
def legendre(l):
    """differentiate the legendre polynomial l times and return result."""
    ex = sp.Symbol("x")
    el = sp.Symbol("l")
    coefficient = 1 / ( 2** el * fact(l) )
    y = (ex**2 - 1)**l
    for i in range(1, l+1):
        y = y.diff(ex)
    return coefficient * y

def assoc_legendre(m, l):
    """differentiate the legendre polynomial m times."""
    ex = sp.Symbol("x")
    m = abs(m)
    original_legendre = legendre(l)
    diff_by_m = original_legendre
    for i in range(1, m+1):
        diff_by_m = diff_by_m.diff(ex)
    function = (1 - ex**2)**(float(m) / 2) * diff_by_m
    # f = sp.lambdify(ex, y)
    # subsitute in cos(x) and return.
    return original_legendre, function

def normalized_angular_solution(m, l):
    """Returns a normalized angular solution for the associated m and l numbers for an orbital."""
    epsilon = -1 ** m if (m >= 0) else 1
    pi_part = ((2 * l) + 1) / (4 * sp.pi)
    separation_factor = fact(1 - abs(m)) / float(fact(1 + abs(m)))
    phi = sp.Symbol("phi")
    original_legendre, diffed_assoc_legendre = assoc_legendre(m, l)
    return original_legendre, diffed_assoc_legendre, epsilon * sp.sqrt(pi_part * separation_factor) * sp.exp(1j*m*phi) * diffed_assoc_legendre
