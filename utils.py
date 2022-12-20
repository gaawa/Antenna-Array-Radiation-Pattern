# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 12:32:41 2022

@author: ThinkWin
"""
import numpy as np

def get_radius_of_equal_distance(d, Npts, ang=2*np.pi):
    """
    Create the radius of the circle curve where "Npts" points is equally distributed
    with distance "d" to the neighboring point where the circle curve has
    the length defined by the angle "ang"

    Parameters
    ----------
    distance : TYPE
        distance between two points.
    Npts : TYPE
        number of points on the curve.
    ang : TYPE, optional
        length of the curve. The default is 2*np.pi.

    Returns
    -------
    radius.

    """
    
    r = d/(2*np.sin( ang/(2*(Npts-1)) ))
    return r
    
def cartesian_to_spherical(x, y, z):
    return np.array([np.sqrt(x**2+y**2+z**2),
                     np.arccos(z/(np.sqrt(x**2+y**2+z**2))),
                     np.sign(y)*np.arccos(x/(x**2+y**2))])

def spherical_to_cartesian(R, theta, phi):
    return np.array([R*np.sin(theta)*np.cos(phi),
                     R*np.sin(theta)*np.sin(phi),
                     R*np.cos(theta)])    