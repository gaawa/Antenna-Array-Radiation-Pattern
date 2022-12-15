# -*- coding: utf-8 -*-
import numpy as np

def cartesian_to_spherical(x, y, z):
    return np.array([np.sqrt(x**2+y**2+z**2),
                     np.arccos(z/(np.sqrt(x**2+y**2+z**2))),
                     np.sign(y)*np.arccos(x/(x**2+y**2))])

def spherical_to_cartesian(R, theta, phi):
    return np.array([R*np.sin(theta)*np.cos(phi),
                     R*np.sin(theta)*np.sin(phi),
                     R*np.cos(theta)])    