import sympy as sp
from sympy.abc import x
from scipy import special
from functions import fact

def c_assoc_laguerre(p, qmp):
    """Generates the associated associated laguerre"""
    return lambda x: int(fact(p + qmp) * special.eval_genlaguerre(qmp, p, x))

def laguerre(q):
    """Generates the laguerre polynomial."""
    coefficient = sp.exp(x)
    arr = []
    y = sp.exp(-x)*x**q
    arr.append(y)
    for i in range(1, q+1):
        y = y.diff(x)
        arr.append(y)
    arr.append(coefficient) # append coefficient into list
    return arr, coefficient * y

def assoc_laguerre(p, qmp):
    main_arr = []
    arr = []
    coefficient = (-1)**p
    q = p + qmp
    laguerre_diffs, diff = laguerre(q)
    main_arr.append(laguerre_diffs)
    arr.append(diff)
    for i in range(1, p+1):
        diff = diff.diff(x)
        arr.append(diff)
    function = coefficient * diff
    main_arr.append(arr)
    return main_arr, function
    