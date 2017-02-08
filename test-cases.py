#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 11:50:13 2017

@author: yaojie
"""

from functions import energy_n, deg_to_rad, rad_to_deg, fact, spherical_to_cartesian, cartesian_to_spherical
import numpy as np

for i in range(1,4):
    print 'n =', i
    print energy_n(i)
print ''

for deg in [90, 180, 270]:
    print 'deg_to_rad({})'.format(deg)
    print deg_to_rad(deg)
print ''

for rad in [3.14, 3.14/2.0, 3.14*3/4]:
    print 'rad_to_deg({})'.format(rad)
    print rad_to_deg(rad)
print ''

for f in [3,5,4,1]:
    print 'fact({})'.format(f)
    print fact(f)
print ''

for r,t,p in [[3,0,np.pi],[3,np.pi/2.0,np.pi/2.0],[3,np.pi,0], [3,60, 25]]:
    print "spherical_to_cartesian({},{},{})".format(r,t,p)
    print spherical_to_cartesian(r,t,p)
print ''

for x,y,z in [[3,0,0],[0,3,0],[0,0,3],[0,-3,0]]:
    print "cartesian_to_spherical({},{},{})".format(x,y,z)
    print cartesian_to_spherical(x,y,z)