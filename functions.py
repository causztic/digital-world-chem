#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 11:48:13 2017

@author: yaojie
"""
import math
import scipy.constants as c
from scipy import special
import numpy as np
import sympy as sp


# week 2

# 1


def energy_n(n):
    ryberg = -c.physical_constants["Rydberg constant times hc in eV"][0]
    return round((ryberg) / (n**2), 5)

# 2


def rad_to_deg(rad):
    return round(float(rad) * (180 / c.pi), 5)


def deg_to_rad(deg):
    return round(float(deg) * (c.pi / 180), 5)

# 3


def spherical_to_cartesian(radius, theta, phi):
    # the theta and phi are swapped, as in the diagram,
    # theta is the elevation instead of phi.
    # Usually, theta is the azimuth.

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(theta)
    return round(x, 5), round(y, 5), round(z, 5)


def cartesian_to_spherical(x, y, z):
    radius = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / radius)
    # use arctan2 instead of arctan to fix values where x is 0.
    phi = np.arctan2(y, x)
    return round(radius, 5), round(theta, 5), round(phi, 5)
########

# week 3

# 4
def fact(f):
    result = 1
    while f > 1:
        result *= f
        f -= 1
    return result