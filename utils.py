# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 12:32:41 2022

@author: ThinkWin
"""
import numpy as np

def get_radius_of_equal_distance(d, Npts, ang=2*np.pi):
    """
    Create the radius of the circle curve along which, "Npts" points are 
    equally distributed with distance "d" to the neighboring point where the 
    circle curve has the length defined by the arc angle "ang"

    Parameters
    ----------
    distance : floating point
        distance between two points.
    Npts : integer
        number of points on the curve.
    ang : floating point, optional
        arc angle of the curve. The default is 2*np.pi (full circle).

    Returns
    -------
    floating point
        radius.

    """
    
    r = d/(2*np.sin( ang/(2*(Npts-1)) ))
    return r
    
def cartesian_to_spherical(x, y, z):
    """
    converts x,y,z values of cartesian coordinate system to a vector containing
    [r,theta,phi] of the spherical coordinate system.

    Parameters
    ----------
    x : floating point / 1D vector
        x-axis values.
    y : floating point / 1D vector
        y-axis value.
    z : floating point / 1D vector
        z-axis value.

    Returns
    -------
    numpy array
        Array containing r, thata and phi of the spherical coordinate system.
        Array dimension depends on the input array dimension.

    """
    return np.array([np.sqrt(x**2+y**2+z**2),
                     np.arccos(z/(np.sqrt(x**2+y**2+z**2))),
                     np.sign(y)*np.arccos(x/(x**2+y**2))])

def spherical_to_cartesian(R, theta, phi):
    """
    converts R,theta,phi values of spherical coordinate system to a vector containing
    [x,y,z] of the cartesian coordinate system.

    Parameters
    ----------
    R : floating point / 1D vector
        Radius values.
    theta : floating point / 1D vector
        Pole angle values.
    phi : floating point / 1D vector
        Azimuth angle values.

    Returns
    -------
    numpy array
        Array containing x, y and z of the cartesian coordinate system.
        Array dimension depends on the input array dimension.

    """
    return np.array([R*np.sin(theta)*np.cos(phi),
                     R*np.sin(theta)*np.sin(phi),
                     R*np.cos(theta)])    