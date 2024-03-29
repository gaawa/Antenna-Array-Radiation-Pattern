# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt


# plot antenna positions
def plot_dots(antennaArray):
    """
    Creates 3D scatter plot of antenna element position

    Parameters
    ----------
    antennaArray : list(AntennaElement)
        List containing antenna elements.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for antennaElement in antennaArray:
        xPos = antennaElement.positionVector[0]
        yPos = antennaElement.positionVector[1]
        zPos = antennaElement.positionVector[2]
        ax.scatter(xPos, yPos, zPos, marker='o', color='b')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()
    
def plot_quiver(antennaArray, ax_size=0.5, arrow_size=1, scale_arrow=1):
    """
    Creatres 3D scatter plot of antenna element position and arrows representing
    antenna element orientation.
    Blue arrow represents vector of 0 pole angle.
    Red arrow represents vector of 0 azimuth angle.

    Parameters
    ----------
    antennaArray : list(AntennaElement)
        List containing antenna elements.
    ax_size : floating point, optional
        Defines axis limits with [-ax_size, ax_size]. The default is 0.5.
    arrow_size : floating point, optional
        Size of arrows for quiver plot. The default is 1.
    scale_arrow : floating point or vector with the same dimension of antennaArray
        scale the arrows belonging to all or to the specific antenna element by the specified amount.

    Returns
    -------
    None.

    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    elementPosList = [antennaElement.positionVector for antennaElement in antennaArray]
    xPosVec, yPosVec, zPosVec = zip(*elementPosList)
    
    thetaVecList = [antennaElement.thetaOrientationVector*arrow_size for antennaElement in antennaArray]
    thetaVecScaled = thetaVecList*scale_arrow
    xDirVecT, yDirVecT, zDirVecT = zip(*thetaVecScaled)
    
    phiVecList = [antennaElement.phiOrientationVector*arrow_size for antennaElement in antennaArray]
    phiVecScaled = phiVecList*scale_arrow
    xDirVecP, yDirVecP, zDirVecP = zip(*phiVecScaled)
    
    # xPosVec = xPosVec+xPosVec
    # yPosVec = yPosVec+yPosVec
    # zPosVec = zPosVec+zPosVec
    # xDirVec = xDirVecT + xDirVecP
    # yDirVec = yDirVecT + yDirVecP
    # zDirVec = zDirVecT + zDirVecP
    
    # ax.quiver(xPosVec, yPosVec, zPosVec, xDirVec, yDirVec, zDirVec)
    
    ax.quiver(xPosVec, yPosVec, zPosVec, xDirVecT, yDirVecT, zDirVecT)
    ax.quiver(xPosVec, yPosVec, zPosVec, xDirVecP, yDirVecP, zDirVecP, color='r')
    ax.scatter(xPosVec, yPosVec, zPosVec, marker='o', color='g')
    
    ax.set_xlim([-ax_size, ax_size])
    ax.set_ylim([-ax_size, ax_size])
    ax.set_zlim([-ax_size, ax_size])
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    plt.show()
    
    