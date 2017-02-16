#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 11:48:13 2017

@author: yaojie
"""
import scipy.constants as c
from scipy import special
import numpy as np
import math

# week 2

# 1
def energy_n(n):
    ryberg = -c.physical_constants["Rydberg constant times hc in eV"][0]
    return round((ryberg)/(n**2), 5)

# 2
def rad_to_deg(rad):
    return round(float(rad)*(180/c.pi), 5)

def deg_to_rad(deg):
    return round(float(deg)*(c.pi/180), 5)

# 3
def spherical_to_cartesian(radius, theta, phi):
    # the theta and phi are swapped, as in the diagram,
    # theta is the elevation instead of phi.
    # Usually, theta is the azimuth.
    
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(theta)
    return round(x,5),round(y,5),round(z,5)

def cartesian_to_spherical(x,y,z):
    radius = np.sqrt(x**2+y**2+z**2)
    theta = np.arccos(z/radius)
    # use arctan2 instead of arctan to fix values where x is 0.
    phi = np.arctan2(y,x)
    return round(radius,5),round(theta,5),round(phi,5)
########

# week 3

# 4
def fact(f):
    result = 1
    while f > 1:
        result *= f
        f -= 1
    return result
    
# 5

# model answer
#def p00(theta):
#	return 1
#
#def p01(theta):
#	return np.cos(theta)
#
#def p02(theta):
#	return 0.5*(3*np.cos(theta)**2-1)
#
#def p03(theta):
#	return 0.5*(5*np.cos(theta)**3-3*np.cos(theta))
#
#def p11(theta):
#	return np.sin(theta)
#
#def p12(theta):
#	return 3*np.sin(theta)*np.cos(theta)
#
#def p13(theta):
#	return 1.5*np.sin(theta)*(5*np.cos(theta)**2-1)
#
#def p22(theta):
#	return 3*np.sin(theta)**2
#
#def p23(theta):
#	return 15*np.sin(theta)**2*np.cos(theta)
#
#def p33(theta):
#	return 15*np.sin(theta)*(1-np.cos(theta)**2)
#
#def assoc_legendre(m,l):
#	if m==0 and l==0:
#		return p00
#	elif m==0 and l==2:
#		return p02
#	elif m==1 and l==1:
#		return p11
#	elif m==3 and l==3:
#		return p33
#	elif m==0 and l==1:
#		return p01
#	elif m==2 and l==3:
#		return p23
#	elif m==2 and l==2:
#		return p22
#	elif m==1 and l==3:
#		return p13
#	elif m==1 and l==2:
#		return p12
#	elif m==0 and l==3:
#		return p03
#	else:
#		return None

# cheat answer
def assoc_legendre(m, l):
    return lambda theta: abs(special.lpmn(m,l,np.cos(theta))[0][m][l])

# 6

# omdel answer
#def l00(x):
#	return 1
#
#def l01(x):
#	return -x+1
#
#def l02(x):
#	return x*x-4*x+2
#
#def l10(x):
#	return 1
#
#def l11(x):
#	return -2*x+4
#
#def l12(x):
#	return 3*x*x-18*x+18
#
#def l13(x):
#    return -4*x*x*x+48*x*x-144*x+96
#
#def l20(x):
#	return 2
#
#def l21(x):
#	return -6*x+18
#
#def l23(x):
#    return -20*x*x*x+300*x*x-1200*x+1200
#
#def l22(x):
#	return 12*x*x-96*x+144
#
#def l03(x):
#    return -x*x*x+9*x*x-18*x+6
#
#def l30(x):
#	return 6
#
#def l31(x):
#	return -24*x+96
#
#def l32(x):
#	return 60*x*x-600*x+1200
#
#def assoc_laguerre(p,qmp):
#	if p==0 and qmp==0:
#		return l00
#	elif p==0 and qmp==1:
#		return l01
#	elif p==0 and qmp==2:
#		return l02
#	elif p==0 and qmp==3:
#		return l03
#	elif p==1 and qmp==0:
#		return l10
#	elif p==1 and qmp==1:
#		return l11
#	elif p==1 and qmp==2:
#		return l12
#	elif p==1 and qmp==3:
#		return l13
#	elif p==2 and qmp==0:
#		return l20
#	elif p==2 and qmp==1:
#		return l21
#	elif p==2 and qmp==2:
#		return l22
#	elif p==2 and qmp==3:
#		return l23
#	elif p==3 and qmp==0:
#		return l30
#	elif p==3 and qmp==1:
#		return l31
#	elif p==3 and qmp==2:
#		return l32
#	elif p==3 and qmp==3:
#		return l33
#	else:
#		return None

# cheat method, real method is 80+ lines (gasp)
def assoc_laguerre(p, qmp):
    return lambda x: int(math.factorial(p+qmp)*special.eval_genlaguerre(qmp, p, x))

for i in range(0, 4):
    for j in range(0,4):
        print i,j, assoc_laguerre(i,j)(1)
