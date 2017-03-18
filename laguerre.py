"""Laguerre Polynomials"""

from math import factorial
import sympy as sp
from sympy.abc import x, a, r, n, l
import scipy.constants as c
from scipy import special

def laguerre(q):
    """Generates and returns the laguerre polynomial with differentials."""
    coefficient = sp.exp(x)
    y = sp.exp(-x) * x**q
    # append the first differential.
    for i in range(1, q + 1):
        y = y.diff(x)
    return coefficient * y

def assoc_laguerre(p, qmp):
    """Generates and returns the associated laguerre polynomial with the previous differentials."""
    coefficient = (-1)**p
    q = p + qmp
    diff = sp.simplify(laguerre(q))
    for i in range(1, p + 1):
        diff = diff.diff(x)
    return coefficient * diff


def normalized_radial_solution(p, qmp):
    """returned the normalized radial solution"""
    coefficient = (2 / (n * a))**3
    coefficient_two = sp.factorial(
        (n - l - 1)) / (2 * n * sp.factorial(n + l) ** 3)
    final_coefficient = sp.sqrt(
        coefficient * coefficient_two) * sp.exp(-r / (n * a)) * ((2 * r) / (n * a))**l

    return final_coefficient * assoc_laguerre(p, qmp)

bohr = c.physical_constants['Bohr radius'][0]

def radial_wave_func(rn, rl, rr):
    p = 2*rl+1
    qmp = rn - rl - 1
    nrs = normalized_radial_solution(p, qmp)
    f = sp.lambdify((n, l, r, a, x), nrs)
    xx = (2 * rr) / (rn * bohr)
    return round(f(rn, rl, rr, bohr, xx) / (bohr ** (-3.0 / 2)), 5)

# print radial_wave_func(1, 0, bohr)
# print radial_wave_func(1, 0, bohr)
# print radial_wave_func(2, 1, bohr)
# print radial_wave_func(2, 1, bohr*2)
# print radial_wave_func(3, 1, bohr*2)