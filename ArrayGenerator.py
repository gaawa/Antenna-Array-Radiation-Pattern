import numpy as np
import coordinate_utils

def generateLinearArray(AntennaElementClass, elementDistance, numElements, xAxisOffset=0):
    antennaArray = []
    
    # thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/4, np.pi/4)
    thetaOrientVec = np.array([0, 0, 1])
    phiOrientVec = np.array([1, 0, 0])
    for n in range(numElements):
        xPos = (elementDistance*n)+xAxisOffset
        positionVector = np.array([xPos, 0, 0])
        antennaArray.append(AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))

    return antennaArray

def generateCircularArray(AntennaElementClass, elementRadius, numElements, poleAng=0, azimuthAngOffset=0):
    angularDistance = 2*np.pi/numElements

    antennaArray = []
    for n in range(numElements):
        elementAzimuth = n*angularDistance + azimuthAngOffset
        xPos = np.cos(elementAzimuth)*elementRadius
        yPos = np.sin(elementAzimuth)*elementRadius
        positionVector = np.array([xPos, yPos, 0])
        thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, poleAng, elementAzimuth)
        phiOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth+np.pi/2)
        antennaArray.append(AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))

    return antennaArray

def generateLinearCircularArray(AntennaElementClass, 
                                elementRadius, numCircularElements,
                                elementLinearDistance, numLinearElements):
    angularDistance = 2*np.pi/numCircularElements
    linearElementPositions = elementLinearDistance*np.array(range(numLinearElements))
    linearArrayLength = np.max(linearElementPositions)
    linearElementCenteredPositions = linearElementPositions - linearArrayLength/2

    antennaArray = []
    for n in range(numCircularElements):
        elementAzimuth = n*angularDistance
        xPos = np.cos(elementAzimuth)*elementRadius
        yPos = np.sin(elementAzimuth)*elementRadius
        for m in range(numLinearElements):
            zPos = linearElementCenteredPositions[m]
            positionVector = np.array([xPos, yPos, zPos])
            # elements are facing outwards
            thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth)
            phiOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth+np.pi/2)
            antennaArray.append(AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))
            
    return antennaArray

def generateBowCircularArray(AntennaElementClass,
                            elementRadius, numCircularElements,
                            elementBowRadius, elementBowAngle, numBowElements):
    
    circAngularDistance = 2*np.pi/numCircularElements

    # generate positions of a arm of bow
    numBowAngSlices = numBowElements-1
    bowAngularDistance = elementBowAngle/numBowAngSlices  # angle between bow elements
    bowElementElevations = -elementBowAngle/2 + np.array(range(numBowElements))*bowAngularDistance
    xPosVec = elementBowRadius * np.cos(bowElementElevations) + (elementRadius - elementBowRadius)
    zPosVec = elementBowRadius * np.sin(bowElementElevations)
    elementPoleAngles = np.pi/2 - bowElementElevations

    antennaArray = []
    for n in range(numCircularElements):
        elementAzimuth = n*circAngularDistance
        rotMatZ = np.array([[np.cos(elementAzimuth), -np.sin(elementAzimuth), 0],
                            [np.sin(elementAzimuth), np.cos(elementAzimuth), 0],
                            [0, 0, 1]])
        for m in range(numBowElements):
            refPosVec = np.array([[xPosVec[m]], [0], [zPosVec[m]]])
            rotPosVec = np.squeeze(np.matmul(rotMatZ, refPosVec))
            thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, 
                                                                     elementPoleAngles[m], 
                                                                     elementAzimuth)
            phiOrientVec = coordinate_utils.spherical_to_cartesian(1, 
                                                                   np.pi/2, 
                                                                   elementAzimuth+np.pi/2)
            antennaArray.append(AntennaElementClass(rotPosVec, thetaOrientVec, phiOrientVec))

    return antennaArray