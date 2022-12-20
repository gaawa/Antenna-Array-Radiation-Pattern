import numpy as np
import coordinate_utils
from utils import get_radius_of_equal_distance, spherical_to_cartesian

class LinearArray():
    def __init__(self, AntennaElementClass, wavelength, elementDistanceFactor=0.5, 
                 numElements=10, xOffset=None):
        self.wavelength = wavelength
        self.elementDistance = wavelength*elementDistanceFactor
        self.numElements = numElements
        if not xOffset:
            self.xOffset = -numElements*self.elementDistance/2
        else:
            self.xOffset = xOffset
        self.AntennaElementClass = AntennaElementClass
        
        self.ax_size = 1.2*numElements*self.elementDistance
        self.arrow_size = 0.5*self.elementDistance
        
        self.arrayElements = self.generateArray()

    def generateArray(self):
        arrayElements = []
        
        # thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/4, np.pi/4)
        thetaOrientVec = np.array([0, 0, 1])
        phiOrientVec = np.array([1, 0, 0])
        for n in range(self.numElements):
            xPos = (self.elementDistance*n)+self.xOffset
            positionVector = np.array([xPos, 0, 0])
            arrayElements.append(self.AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))
    
        return arrayElements
    
class CircularArray():
    def __init__(self, AntennaElementClass, wavelength, elementDisctanceFactor=0.5, 
                 numElements=10, circularAng=2*np.pi, 
                 elementPoleAng=np.pi/2, circularAzimuthOffset=0):
        self.wavelength = wavelength
        self.elementDistance = wavelength*elementDisctanceFactor
        self.numElements = numElements
        self.circularAng = circularAng
        self.AntennaElementClass = AntennaElementClass
        self.elementRadius = get_radius_of_equal_distance(self.elementDistance, 
                                                          numElements,
                                                          ang=circularAng)
        self.elementPoleAng = elementPoleAng
        self.circularAzimuthOffset = circularAzimuthOffset
        
        self.ax_size = 1.2*2*self.elementRadius
        self.arrow_size = 0.5*(2*np.pi*self.elementRadius/numElements)
        
        self.arrayElements = self.generateArray()


    def generateArray(self):
        angularDistance = 2*np.pi/self.numElements
    
        arrayElements = []
        for n in range(self.numElements):
            elementAzimuth = n*angularDistance + self.circularAzimuthOffset
            xPos = np.cos(elementAzimuth)*self.elementRadius
            yPos = np.sin(elementAzimuth)*self.elementRadius
            positionVector = np.array([xPos, yPos, 0])
            thetaOrientVec = spherical_to_cartesian(1, self.elementPoleAng, elementAzimuth)
            phiOrientVec = spherical_to_cartesian(1, np.pi/2, elementAzimuth+np.pi/2)
            arrayElements.append(self.AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))
    
        return arrayElements
    
class LinearCircularArray():
    def __init__(self, AntennaElementClass, wavelength,
                 linearDistanceFactor=0.5, circularDistanceFactor=0.5,
                 numLinearElements=10, numCircularElements=10,
                 circularAng=2*np.pi):
        self.AntennaElementClass = AntennaElementClass
        self.wavelength = wavelength
        self.linearElementDistance = wavelength*linearDistanceFactor
        self.circularElementDistance = wavelength*circularDistanceFactor
        self.numLinearElements = numLinearElements
        self.numCircularElements = numCircularElements
        self.circularAng = circularAng
        self.circularRadius = get_radius_of_equal_distance(self.circularElementDistance, 
                                                           self.numCircularElements,
                                                           ang=circularAng)
        
        self.ax_size = 1.2 * max(self.numLinearElements*self.linearElementDistance,
                        2*self.circularRadius)
        self.arrow_size = 0.5 * min(self.linearElementDistance,
                               2*np.pi*self.circularRadius/self.numCircularElements)
        
        self.arrayElements = self.generateArray()
        
    def generateArray(self):
        angularDistance = 2*np.pi/self.numCircularElements
        linearElementPositions = self.linearElementDistance*np.array(range(self.numLinearElements))
        linearArrayLength = np.max(linearElementPositions)
        linearElementCenteredPositions = linearElementPositions - linearArrayLength/2
    
        arrayElements = []
        for n in range(self.numCircularElements):
            elementAzimuth = n*angularDistance
            xPos = np.cos(elementAzimuth)*self.circularRadius
            yPos = np.sin(elementAzimuth)*self.circularRadius
            for m in range(self.numLinearElements):
                zPos = linearElementCenteredPositions[m]
                positionVector = np.array([xPos, yPos, zPos])
                # elements are facing outwards
                thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth)
                phiOrientVec = coordinate_utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth+np.pi/2)
                arrayElements.append(self.AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))
                
        return arrayElements  
    
class BowCircularArray():
    def __init__(self, AntennaElementClass, wavelength,
                 bowDistanceFactor=0.5, circularDistanceFactor=0.5,
                 numBowElements=10, numCircularElements=10,
                 bowAng=np.pi/2, circularAng=2*np.pi):
        self.AntennaElementClass = AntennaElementClass
        self.wavelength = wavelength
        self.bowElementDistance = wavelength*bowDistanceFactor
        self.circularElementDistance = wavelength*circularDistanceFactor
        self.numBowElements = numBowElements
        self.numCircularElements = numCircularElements
        self.bowAng = bowAng
        self.bowRadius = get_radius_of_equal_distance(self.bowElementDistance, 
                                                      numBowElements,
                                                      ang=self.bowAng)
        self.circularAng = circularAng
        if numCircularElements == 1:
            self.circularRadius = 0.1
        else:
            self.circularRadius = get_radius_of_equal_distance(self.circularElementDistance, 
                                                           numCircularElements,
                                                           ang=self.circularAng)
        
        self.ax_size = 1.2 * 2*self.circularRadius
        arrow_size_circ = 0.5 * 2*np.pi*self.circularRadius/numCircularElements
        arrow_size_bow = 0.5 * 2*np.pi*self.bowRadius/numBowElements
        self.arrow_size = min(arrow_size_circ, arrow_size_bow)
        
        self.arrayElements = self.generateArray()
        
    def generateArray(self):
        circAngularDistance = 2*np.pi/self.numCircularElements
    
        # generate positions of a arm of bow
        numBowAngSlices = self.numBowElements-1
        bowAngularDistance = self.bowAng/numBowAngSlices  # angle between bow elements
        bowElementElevations = -self.bowAng/2 + np.array(range(self.numBowElements))*bowAngularDistance
        xPosVec = self.bowRadius * np.cos(bowElementElevations) + (self.circularRadius - self.bowRadius)
        zPosVec = self.bowRadius * np.sin(bowElementElevations)
        elementPoleAngles = np.pi/2 - bowElementElevations
    
        arrayElements = []
        for n in range(self.numCircularElements):
            elementAzimuth = n*circAngularDistance
            rotMatZ = np.array([[np.cos(elementAzimuth), -np.sin(elementAzimuth), 0],
                                [np.sin(elementAzimuth), np.cos(elementAzimuth), 0],
                                [0, 0, 1]])
            for m in range(self.numBowElements):
                refPosVec = np.array([[xPosVec[m]], [0], [zPosVec[m]]])
                rotPosVec = np.squeeze(np.matmul(rotMatZ, refPosVec))
                thetaOrientVec = coordinate_utils.spherical_to_cartesian(1, 
                                                                         elementPoleAngles[m], 
                                                                         elementAzimuth)
                phiOrientVec = coordinate_utils.spherical_to_cartesian(1, 
                                                                       np.pi/2, 
                                                                       elementAzimuth+np.pi/2)
                arrayElements.append(self.AntennaElementClass(rotPosVec, thetaOrientVec, phiOrientVec))
    
        return arrayElements