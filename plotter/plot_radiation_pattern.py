# -*- coding: utf-8 -*-

# https://stackoverflow.com/questions/54822873/python-plotting-antenna-radiation-pattern;

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import utils

def plot2d(arrayPattern, poleAngles, azimuthAngles):
    """
    2D surf plot of radiation pattern.

    Parameters
    ----------
    arrayPattern : 2D matrix
        Matrix containing signal amplification correspondig to poleAngles and
        azimuthAngles.
    poleAngles : 1D matrix
        pole angle values.
    azimuthAngles : 1D matrix
        azimuth angle values.

    Returns
    -------
    None.

    """
    arrayPatternLog = 10*np.log10(np.abs((arrayPattern)**2))
    arrayPatternLog[arrayPatternLog<-30] = np.nan
    
    # plot radiation pattern
    iMesh, jMesh = np.meshgrid(poleAngles, azimuthAngles, indexing='ij')
    iMesh = iMesh/(2*np.pi)*360
    jMesh = jMesh/(2*np.pi)*360
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(jMesh, iMesh, arrayPatternLog, cmap=cm.jet)
    ax.set_xlabel('phi')
    ax.set_ylabel('theta')
    ax.set_zlim(np.nanmax(arrayPatternLog)-40, np.nanmax(arrayPatternLog))
    plt.show()

def plot3d(arrayPattern, poleAngles, azimuthAngles, logMode=True, logRange = 40):
    """
    3D surf plot of radiation pattern.

    Parameters
    ----------
    arrayPattern : 2D matrix
        Matrix containing signal amplification correspondig to poleAngles and
        azimuthAngles.
    poleAngles : 1D matrix
        pole angle values.
    azimuthAngles : 1D matrix
        azimuth angle values.
    logMode : boolean, optional
        Values are converted with 20*log10(abs(arrayPattern)) first if True. 
        The default is True.
    logRange : floating point, optional
        The value range from maximum value downto max-logRange that are plotted. 
        Smaller values are ignored.
        The default is 40.

    Returns
    -------
    None.

    """
    if logMode:
        # set values below the limit of max-logRange to the limit
        arrayPatternLog = 10*np.log10(np.abs((arrayPattern)**2))
        correction = -arrayPatternLog.max()+logRange
        arrayPatternLog = arrayPatternLog + correction
        arrayPatternLog[arrayPatternLog<0] = 0
        
        iMesh, jMesh = np.meshgrid(poleAngles, azimuthAngles, indexing='ij')
        arrayPatternCoord = utils.spherical_to_cartesian(arrayPatternLog, 
                                                                    iMesh, jMesh)
        R = arrayPatternLog
        
    else:
        # TODO: scalar plot
        pass
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.grid(True)
    ax.axis('on')
    
    mycol = cm.jet(R/np.max(R))
    surf = ax.plot_surface(arrayPatternCoord[0,:], 
                           arrayPatternCoord[1,:], 
                           arrayPatternCoord[2,:], 
                           rstride=1, cstride=1, facecolors=mycol, 
                           linewidth=0.5, antialiased=True, shade=False, 
                           alpha=0.5, zorder = 0.5)
    
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(R-correction)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    Rmax = np.max(R)
    ax.set_xlim([-Rmax, Rmax])
    ax.set_ylim([-Rmax, Rmax])
    ax.set_zlim([-Rmax, Rmax])
    
    fig.colorbar(m, shrink=0.8)
    ax.view_init(azim=300, elev=30)
    plt.show()