# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt


# plot antenna positions
def plotDots(antennaArray):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for antennaElement in antennaArray:
        xPos = antennaElement.positionVector[0]
        yPos = antennaElement.positionVector[1]
        zPos = antennaElement.positionVector[2]
        ax.scatter(xPos, yPos, zPos, marker='o', color='b')
    plt.show()
    
def plotQuiver(antennaArray, ax_size=0.5, arrow_size=1):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    elementPosList = [antennaElement.positionVector for antennaElement in antennaArray]
    xPosVec, yPosVec, zPosVec = zip(*elementPosList)
    
    thetaVecList = [antennaElement.thetaOrientationVector*arrow_size for antennaElement in antennaArray]
    xDirVecT, yDirVecT, zDirVecT = zip(*thetaVecList)
    
    phiVecList = [antennaElement.phiOrientationVector*arrow_size for antennaElement in antennaArray]
    xDirVecP, yDirVecP, zDirVecP = zip(*phiVecList)
    
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
    
    