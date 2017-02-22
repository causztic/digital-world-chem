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
    l = sp.exp(-x)*x**q
    for i in range(1, q+1):
        l = l.diff(x)
        