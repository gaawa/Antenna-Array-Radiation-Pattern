import numpy as np
import utils
from utils import get_radius_of_equal_distance, spherical_to_cartesian
from abc import ABC


class AntennaArray(ABC):
    """
    Abstract base class for all antenna array classes
    """
    def __init__(self, arrayElements, ax_size, arrow_size):
        """
        This ABC requires the list of antenna elements, axis size for antenna 
        element plot by [-ax_size, ax_size] and the arrow size for overlaying 
        quiver plot.

        Parameters
        ----------
        AntennaElementClass : AntennaElement
            The AntennaElement class from which to instantiate each antenna element from.
        arrayElements : list(AntennaElement)
            list of antenna element objects belonging to this antenna array.
        ax_size : floating point
            axis size for antenna element plot by [-ax_size, ax_size].
        arrow_size : floating point
            arrow
            size for overlaying quiver plot.

        Returns
        -------
        None.

        """
        self.arrayElements = arrayElements
        self.ax_size = ax_size
        self.arrow_size = arrow_size

class LinearArray(AntennaArray):
    """
    Class that creates and stores attributes of a linear array formation
    """
    def __init__(self, AntennaElementClass, wavelength, elementDistanceFactor=0.5, 
                 numElements=10, xOffset=None):
        """
        Create linear antenna array.

        Parameters
        ----------
        AntennaElementClass : AntennaElement
            The AntennaElement class from which to instantiate each antenna element from.
        wavelength : floating point
            lambda of the system.
        elementDistanceFactor : floating poing, optional
            Defines the distance of the antenna element.
            elementDistance = wavelength*elementDistanceFactor.
            The default is 0.5.
        numElements : integer, optional
            number of elements for the linear array. The default is 10.
        xOffset : floating point, optional
            x-position offset of the linear array position. The default is None.

        Returns
        -------
        None.

        """
        self.AntennaElementClass = AntennaElementClass
        self.wavelength = wavelength
        self.elementDistance = wavelength*elementDistanceFactor
        self.numElements = numElements
        if not xOffset:
            self.xOffset = -numElements*self.elementDistance/2
        else:
            self.xOffset = xOffset
            
        arrayElements = self.generateArray()
        
        ax_size = 1.2*numElements*self.elementDistance
        arrow_size = 0.5*self.elementDistance
        
        super().__init__(arrayElements, ax_size, arrow_size)

    def generateArray(self):
        arrayElements = []
        
        # thetaOrientVec = utils.spherical_to_cartesian(1, np.pi/4, np.pi/4)
        thetaOrientVec = np.array([0, 0, 1])
        phiOrientVec = np.array([1, 0, 0])
        for n in range(self.numElements):
            xPos = (self.elementDistance*n)+self.xOffset
            positionVector = np.array([xPos, 0, 0])
            arrayElements.append(self.AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))
    
        return arrayElements
    
class CircularArray(AntennaArray):
    """
    Class that creates and stores attributes of a circular array formation
    """
    def __init__(self, AntennaElementClass, wavelength, elementDisctanceFactor=0.5, 
                 numElements=10, circularAng=2*np.pi, 
                 elementPoleAng=np.pi/2, circularAzimuthOffset=0):
        """
        Create circular array.

        Parameters
        ----------
        AntennaElementClass : AntennaElement
            The AntennaElement class from which to instantiate each antenna element from.
        wavelength : floating point
            lambda of the system.
        elementDisctanceFactor : floating point, optional
            Defines the distance of the antenna element.
            elementDistance = wavelength*elementDistanceFactor.
            The default is 0.5.
        numElements : integer, optional
            number of elements for the circular array. The default is 10.
        circularAng : floating poing, optional
            Angle of the arc. The default is 2*np.pi (whole circle).
        elementPoleAng : floating poing, optional
            Pole angle orientation of the antenna elements. 
            The default is np.pi/2 (pointing outwards).
        circularAzimuthOffset : floating poing, optional
            Azimuth angle offset where the first antenna element shall be positioned. 
            The default is 0.

        Returns
        -------
        None.

        """
        self.AntennaElementClass = AntennaElementClass
        self.wavelength = wavelength
        self.elementDistance = wavelength*elementDisctanceFactor
        self.numElements = numElements
        self.circularAng = circularAng
        self.elementRadius = get_radius_of_equal_distance(self.elementDistance, 
                                                          numElements,
                                                          ang=circularAng)
        self.elementPoleAng = elementPoleAng
        self.circularAzimuthOffset = circularAzimuthOffset
        
        ax_size = 1.2*2*self.elementRadius
        arrow_size = 0.5*(2*np.pi*self.elementRadius/numElements)
        
        arrayElements = self.generateArray()
        
        super().__init__(arrayElements, ax_size, arrow_size)


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
    
class LinearCircularArray(AntennaArray):
    """
    Class that creates and stores attributes of a cylindrical array formation
    """
    def __init__(self, AntennaElementClass, wavelength,
                 linearDistanceFactor=0.5, circularDistanceFactor=0.5,
                 numLinearElements=4, numCircularElements=10,
                 circularAng=2*np.pi):
        """
        Create cylindrical array.
        The cylindrical array consists of multiple linear arrays arranged in a circle.

        Parameters
        ----------
        AntennaElementClass : AntennaElement
            The AntennaElement class from which to instantiate each antenna element from.
        wavelength : floating poing
            lambda of the system.
        linearDistanceFactor : floating point, optional
            Defines the distance of the linear array antenna element.
            linearElementDistance = wavelength*linearDistanceFactor.
            The default is 0.5.
        circularDistanceFactor : floating point, optional
            Defines the distance of the circular array antenna element.
            circularElementDistance = wavelength*circularDistanceFactor.
            The default is 0.5.
        numLinearElements : int, optional
            number of elements for the linear array. The default is 10.
        numCircularElements : TYPE, optional
            number of elements for the circular array. The default is 10.
        circularAng : TYPE, optional
            Angle of the circle arc. The default is 2*np.pi (whole circle).

        Returns
        -------
        None.

        """
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
        
        ax_size = 1.2 * max(self.numLinearElements*self.linearElementDistance,
                            2*self.circularRadius)
        arrow_size = 0.5 * min(self.linearElementDistance,
                               2*np.pi*self.circularRadius/self.numCircularElements)
        
        arrayElements = self.generateArray()
        
        super().__init__(arrayElements, ax_size, arrow_size)
        
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
                thetaOrientVec = utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth)
                phiOrientVec = utils.spherical_to_cartesian(1, np.pi/2, elementAzimuth+np.pi/2)
                arrayElements.append(self.AntennaElementClass(positionVector, thetaOrientVec, phiOrientVec))
                
        return arrayElements  
    
class BowCircularArray(AntennaArray):
    """
    Clas that creates and stores attributes of a bow circular array formation
    """
    def __init__(self, AntennaElementClass, wavelength,
                 bowDistanceFactor=0.5, circularDistanceFactor=0.5,
                 numBowElements=4, numCircularElements=10,
                 bowAng=np.pi/2, circularAng=2*np.pi):
        """
        Create bow-circular array.
        Bow-circular array consists of multiple bow arrays (antennas placed along
        an arc) arranged in a circle.

        Parameters
        ----------
        AntennaElementClass : AntennaElement
            The AntennaElement class from which to instantiate each antenna element from.
        wavelength : floating point
            lambda of the system.
        bowDistanceFactor : floating point, optional
            Defines the distance of the bow array antenna elements.
            bowElementDistance = wavelength*bowDistanceFactor.
            The default is 0.5.
        circularDistanceFactor : floating point, optional
            Defines the distance of the circular array antenna element.
            circularElementDistance = wavelength*circularDistanceFactor.
            The default is 0.5.
        numBowElements : integer, optional
            number of elements for the bow array. The default is 10.
        numCircularElements : integer, optional
            number of elements for the circular array. The default is 10.
        bowAng : TYPE, optional
            Angle of the bow arc. The default is np.pi/2 (half circle).
        circularAng : TYPE, optional
            Angle of the circle arc. The default is 2*np.pi (whole circle).

        Returns
        -------
        None.

        """
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
        
        ax_size = 1.2 * 2*self.circularRadius
        arrow_size_circ = 0.5 * 2*np.pi*self.circularRadius/numCircularElements
        arrow_size_bow = 0.5 * 2*np.pi*self.bowRadius/numBowElements
        arrow_size = min(arrow_size_circ, arrow_size_bow)
        
        arrayElements = self.generateArray()
        
        super().__init__(arrayElements, ax_size, arrow_size)
        
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
                thetaOrientVec = utils.spherical_to_cartesian(1, 
                                                                         elementPoleAngles[m], 
                                                                         elementAzimuth)
                phiOrientVec = utils.spherical_to_cartesian(1, 
                                                                       np.pi/2, 
                                                                       elementAzimuth+np.pi/2)
                arrayElements.append(self.AntennaElementClass(rotPosVec, thetaOrientVec, phiOrientVec))
    
        return arrayElements