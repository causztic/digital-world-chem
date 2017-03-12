import numpy
import scipy.constants as c
import sympy as sp
import math
from math import factorial
from sympy.abc import x, a, r, n, l

def cartesian_to_spherical(x, y, z):
    radius = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / radius)
    # use arctan2 instead of arctan to fix values where x is 0.
    phi = np.arctan2(y, x)
    return round(radius, 5), round(theta, 5), round(phi, 5)

def laguerre(q):
    """Generates and returns the laguerre polynomial with differentials."""
    coefficient = sp.exp(x)
    arr = []
    y = sp.exp(-x) * x**q
    # append the first differential.
    arr.append(y)
    for i in range(1, q + 1):
        y = y.diff(x)
        # for every differential, append it to arr.
        arr.append(y)
    arr.append(coefficient)  # append coefficient into list
    return arr, coefficient * y

def assoc_laguerre(p, qmp):
    """Generates and returns the associated laguerre polynomial with the previous differentials."""
    main_arr = []
    arr = []
    coefficient = (-1)**p
    q = p + qmp
    laguerre_diffs, diff = laguerre(q)
    diff = sp.simplify(diff)
    main_arr.append(laguerre_diffs)
    arr.append(diff)
    for i in range(1, p + 1):
        diff = diff.diff(x)
        arr.append(diff)
    function = coefficient * diff
    main_arr.append(arr)
    return main_arr, function


def normalized_radial_solution(p, qmp):
    """returned the normalized radial solution"""
    diffs, diff_assoc_laguerre = assoc_laguerre(p, qmp)
    coefficient = (2 / (n * a))**3
    coefficient_two = sp.factorial(
        (n - l - 1)) / (2 * n * sp.factorial(n + l) ** 3)
    final_coefficient = sp.sqrt(
        coefficient * coefficient_two) * sp.exp(-r / (n * a)) * ((2 * r) / (n * a))**l

    return diffs, diff_assoc_laguerre, final_coefficient * diff_assoc_laguerre

bohr = c.physical_constants['Bohr radius'][0]

def radial_wave_func(rn, rl, rr):
    p = 2*rl+1
    qmp = rn - rl - 1
    nrs = normalized_radial_solution(p, qmp)[2]
    f = sp.lambdify((n, l, r, a, x), nrs)
    xx = (2 * rr) / (rn * bohr)
    return round(f(rn, rl, rr, bohr, xx) / (bohr ** (-3.0 / 2)), 5)

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

def hydrogen_wave_func(n, l, m, roa, nx, ny, nz):
    roa = float(roa)
    xx, yy, zz = numpy.meshgrid(numpy.linspace(-roa, roa, nx),
                                numpy.linspace(-roa, roa, ny), numpy.linspace(-roa, roa, nz))

    vector_r = numpy.vectorize(cartesian_to_spherical)
    v_angularwf = numpy.vectorize(angular_wave_func)
    v_radialwf = numpy.vectorize(radial_wave_func)

    bohr_radius = c.physical_constants['Bohr radius'][0]
    r, theta, phi = vector_r(xx, yy, zz)
    r = r * bohr_radius

    radial_part = v_radialwf(n, l, r)
    angular_part = v_angularwf(m, l, theta, phi)

    rwf = (radial_part * angular_part) ** 2
    mag = numpy.absolute(rwf)

    return numpy.round(xx,5), numpy.round(yy,5), numpy.round(zz,5), numpy.round(mag,5)


