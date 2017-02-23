#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 11:50:13 2017

@author: yaojie
"""

# from functions import energy_n, deg_to_rad, rad_to_deg, fact, spherical_to_cartesian, cartesian_to_spherical
from laguerre import c_assoc_laguerre, assoc_laguerre, normalized_radial_solution
from legendre import assoc_legendre, c_assoc_legendre, normalized_angular_solution
from model_answers import assoc_laguerre as m_assoc_laguerre, assoc_legendre as m_assoc_legendre
import sympy as sp

RED_START = '\33[31m'
END = '\033[0m'
print "Associated Legendre: "

s = ""
for l in range(0, 4):
    for m in range(-l, l + 1):
        theta = 1
        mal = m_assoc_legendre(m, l)
        cal = c_assoc_legendre(m, l)
        if mal != None and cal != None:
            s += "%sm\t l\t Θ\tmodel\t cheat\t actual\t yay?%s\t\n" % (RED_START, END)
            mal = mal(theta)
            cal = cal(theta)
            aal = assoc_legendre(m, l)[1]
            [l_diffs, m_diffs], assoc_func, normal_func = normalized_angular_solution(m, l)
            actual = sp.lambdify(
                (sp.Symbol('x'), sp.Symbol("l")), aal)(sp.cos(theta), l)
            check = 'boo'
            if round(mal, 9) == round(actual, 9):
                check = 'woo!'
            s += "%s\t %d\t %s\t" % (m, l, theta)
            s += "%.3f\t %.3f\t %.3f\t %s\t\n" % (mal, cal, actual, check)
            s += "legendre differentiation: \n"
            for idx, i in enumerate(l_diffs):
                if idx == len(l_diffs)-1:
                    s += "Adding in coefficient => %s\n" % i
                else:
                    # legendre polynomial diffs with l.
                    s += "differentiation l = #%d %s\n" % (idx, i)

            for idx, i in enumerate(m_diffs):
                # associated legendre polynomial diffs with m.
                s += "differentiation m = #%d %s\n" % (idx, i)

            s += "associated legendre function: \t%s\n" % assoc_func
            s += "normalized angular solution: \t%s\n\n" % normal_func
print s
s = ""

print "Associated Laguerre: "
for n in range(1, 5):
    for l in range(0, n):
        p = 2 * l + 1
        qmp = n - l - 1
        s += "\nn = %i, l = %i\n" % (n, l)
        x = 1
        mal = m_assoc_laguerre(p, qmp)
        cal = c_assoc_laguerre(p, qmp)

        s += "\n%sp qmp x\tmodel\tcheat\t actual\t yay?%s\t\n" % (RED_START, END)
        cal = cal(x)
        aal = assoc_laguerre(p, qmp)[1]
        actual = sp.lambdify(sp.Symbol('x'), aal)(x)
        [q_diffs, p_diffs], assoc_func, normal_func = normalized_radial_solution(p, qmp)
        if mal != None:
            mal = mal(x)
        else:
            mal = cal

        check = 'boo'
        if abs(mal - actual) <= 10**-9:
            check = 'woo!'
        s += "%s %d   %s\t" % (p, qmp, x)
        s += "%.1f\t%.1f\t%.1f\t %s\t\n" % (mal, cal, actual, check)
        s += "laguerre differentiation: \n"
        for idx, i in enumerate(q_diffs):
            if idx == len(q_diffs)-1:
                s += "Adding in coefficient => %s\n" % i
            else:
                # laguerre polynomial differentation with q
                s += "differentiation q = #%d %s\n" % (idx, i)

        for idx, i in enumerate(p_diffs):
            # associated laguerre differentation with qmp
            s += "differentiation p = #%d %s\n" % (idx, i)
        
        s += "associated laguerre function: \t%s\n" % assoc_func
        s += "normalized radial solution: \t%s\n\n" % normal_func
print s
# for i in range(1,4):
#     print 'n =', i
#     print energy_n(i)
# print ''

# for deg in [90, 180, 270]:
#     print 'deg_to_rad({})'.format(deg)
#     print deg_to_rad(deg)
# print ''

# for rad in [3.14, 3.14/2.0, 3.14*3/4]:
#     print 'rad_to_deg({})'.format(rad)
#     print rad_to_deg(rad)
# print ''

# for f in [3,5,4,1]:
#     print 'fact({})'.format(f)
#     print fact(f)
# print ''

# for r,t,p in [[3,0,np.pi],[3,np.pi/2.0,np.pi/2.0],[3,np.pi,0], [3,60, 25]]:
#     print "spherical_to_cartesian({},{},{})".format(r,t,p)
#     print spherical_to_cartesian(r,t,p)
# print ''

# for x,y,z in [[3,0,0],[0,3,0],[0,0,3],[0,-3,0]]:
#     print "cartesian_to_spherical({},{},{})".format(x,y,z)
#     print cartesian_to_spherical(x,y,z)
