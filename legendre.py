"""
Contains the Legendre Polynomial calculation for Schrodinger's Equation.
The return values of each function will return a symbolic equation which can be
lambdified with Sympy to generate a solution. The differential equations will
also be removed for posterity.
"""

from math import factorial
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
    coefficient = 1 / (2 ** el * factorial(l))
    arr = []
    y = (ex**2 - 1)**l
    arr.append(y)  # add non-differentiated equation
    for i in range(1, l + 1):
        y = y.diff(ex)
        # for every differentiation, append new equation.
        arr.append(y)
    arr.append(coefficient)
    return arr, coefficient * y


def assoc_legendre(m, l):
    """differentiate the legendre polynomial m times."""
    main_arr = []
    arr = []
    ex = sp.Symbol("x")
    m = abs(m)
    legendre_diffs, diff_by_m = legendre(l)
    main_arr.append(legendre_diffs)
    arr.append(diff_by_m)
    for i in range(1, m + 1):
        diff_by_m = diff_by_m.diff(ex)
        # append additional differentiations by m per m differentiaion.
        arr.append(diff_by_m)
    function = (1 - ex**2)**(float(m) / 2) * diff_by_m
    # append assoc differenations by m, if any.
    main_arr.append(arr)
    return main_arr, function


def normalized_angular_solution(m, l):
    """Returns a normalized angular solution for the associated m and l numbers for an orbital."""
    epsilon = -1 ** m if (m >= 0) else 1
    pi_part = ((2 * l) + 1) / (4 * sp.pi)
    separation_factor = factorial(l - abs(m)) / float(factorial(l + abs(m)))
    phi = sp.Symbol("phi")
    diffs, diffed_assoc_legendre = assoc_legendre(m, l)
    funct =  epsilon * sp.sqrt(pi_part * separation_factor) * sp.exp(1j * sp.im(m) * sp.im(phi)) * diffed_assoc_legendre
    return diffs, diffed_assoc_legendre, funct

def angular_wave_func(m, l, theta, phi):
    nas = normalized_angular_solution(m, l)[2]
    f = sp.lambdify((sp.Symbol("x"), sp.Symbol("phi"), sp.Symbol("l")), nas)
    return complex(round(-1*f(np.cos(theta), phi, l),5))

# print angular_wave_func(0,0,0,0)
# print angular_wave_func(0,1,np.pi,0)
# print angular_wave_func(1,1,np.pi/2,np.pi)
# print angular_wave_func(0,2,np.pi,0)