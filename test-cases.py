#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 11:50:13 2017

@author: yaojie
"""

# from functions import energy_n, deg_to_rad, rad_to_deg, fact, spherical_to_cartesian, cartesian_to_spherical
from functions import c_assoc_laguerre
from legendre import assoc_legendre, c_assoc_legendre, normalized_angular_solution
from model_answers import assoc_laguerre as m_assoc_laguerre, assoc_legendre as m_assoc_legendre
import numpy as np
import sympy as sp
RED_START = '\33[31m'
END = '\033[0m'
print "Associated Legendre: "
s = ""
for l in range(0, 4):
    for m in range(-l, l + 1):
        for theta in [0, 1]:
            mal = m_assoc_legendre(m, l)
            cal = c_assoc_legendre(m, l)
            if mal != None and cal != None:
                s += "%sm\t l\t Θ\tmodel\t cheat\t actual\t yay?%s\t\n" % (RED_START, END)
                mal = mal(theta)
                cal = cal(theta)
                aal = assoc_legendre(m, l)[1]
                nas = normalized_angular_solution(m, l)
                actual = sp.lambdify(
                    (sp.Symbol('x'), sp.Symbol("l")), aal)(sp.cos(theta), l)
                check = 'boo'
                if round(mal, 9) == round(actual, 9):
                    check = 'woo!'
                s += "%s\t %d\t %s\t" % (m, l, theta)
                s += "%.3f\t %.3f\t %.3f\t %s\t\n" % (mal, cal, actual, check)
                s += "legendre differentiation: \t\t\t\t%s\n" % nas[0]
                s += "associated legendre differentiation: \t%s\n" % nas[1]
                s += "normalized angular solution: \t\t\t%s\n\n" % nas[2]
print s

print "Associated Laguerre: "
for p in range(0, 4):
    for qmp in range(0, 4):
        for x in [0, 1]:
            mal = m_assoc_laguerre(p, qmp)
            cal = c_assoc_laguerre(p, qmp)
            if mal != None and cal != None:
                print "p = %s\t qmp = %s\t x = %s\t %s\t %s" % (p, qmp, x, mal(x), cal(x))

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
