import sympy as sp
from sympy.abc import x,a,r,n,l
from scipy import special
from functions import fact

def c_assoc_laguerre(p, qmp):
    """Generates the cheat associated laguerre"""
    return lambda x: int(fact(p + qmp) * special.eval_genlaguerre(qmp, p, x))

def laguerre(q):
    """Generates and returns the laguerre polynomial with differentials."""
    coefficient = sp.exp(x)
    arr = []
    y = sp.exp(-x)*x**q
    # append the first differential.
    arr.append(y)
    for i in range(1, q+1):
        y = y.diff(x)
        # for every differential, append it to arr.
        arr.append(y)
    arr.append(coefficient) # append coefficient into list
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
    for i in range(1, p+1):
        diff = diff.diff(x)
        arr.append(diff)
    function = coefficient * diff
    main_arr.append(arr)
    return main_arr, function

def normalized_radial_solution(p, qmp):

    diffs, diff_assoc_laguerre = assoc_laguerre(p, qmp)
    coefficient = (2/(n * a))**3
    coefficient_two = sp.factorial((n-l-1)) / (2 * n * sp.factorial(n+l) ** 3)
    final_coefficient = sp.sqrt(coefficient * coefficient_two) * sp.exp(-r/(n*a)) * ((2 * r) / (n * a))**l

    return diffs, diff_assoc_laguerre, final_coefficient * diff_assoc_laguerre
