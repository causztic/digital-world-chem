#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 11:48:13 2017

@author: yaojie
"""
import scipy.constants as c
from scipy import special
import numpy as np

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
def assoc_laguerre(q, p, x):
    return special.assoc_laguerre(x, p,q-p)
    
print assoc_laguerre(0,0,1)
print assoc_laguerre(1,1,1)
print assoc_laguerre(2,2,1)
print assoc_laguerre(2,2,0)
