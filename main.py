import matplotlib.pyplot as plt
from matplotlib import cm
import ArrayGenerator
from OmnidirectionalAntennaElement import OmnidirectionalAntennaElement
import scipy.constants as spc
import numpy as np
from enum import Enum, auto

class ArrayFormation(Enum):
    Linear = auto()
    Circular = auto()
    LinearCircular = auto()
    BowCircular = auto()

# spawn antenna elements in linear array
freq = 3e9
wavelength = spc.speed_of_light/freq
wavenumber = 2*np.pi/wavelength
elementDistance = wavelength/2
arrayFormation = ArrayFormation.LinearCircular

if arrayFormation is ArrayFormation.Linear:
    numElements = 10
    antennaArray = ArrayGenerator.generateLinearArray(OmnidirectionalAntennaElement,
                                                      elementDistance, numElements)

elif arrayFormation is ArrayFormation.Circular:
    numElements = 10
    elementRadius = wavelength/(4*np.sin(np.pi/numElements))
    antennaArray = ArrayGenerator.generateCircularArray(OmnidirectionalAntennaElement,
                                                        elementRadius, numElements)

elif arrayFormation is ArrayFormation.LinearCircular:
    numCircularElements = 10
    elementRadius = wavelength/(4*np.sin(np.pi/numCircularElements))
    numLinearElements = 4
    elementLinearDistance = wavelength/2
    antennaArray = ArrayGenerator.generateLinearCircularArray(OmnidirectionalAntennaElement,
                                                              elementRadius, numCircularElements,
                                                              elementLinearDistance, numLinearElements)

elif arrayFormation is ArrayFormation.BowCircular:
    numCircularElements = 10
    elementRadius = wavelength/(4*np.sin(np.pi/numCircularElements))
    numBowElements = 4
    elementBowAngle = np.pi/2
    elementBowRadius = elementDistance/(2*np.sin(elementBowAngle/(2*numBowElements)))
    antennaArray = ArrayGenerator.generateBowCircularArray(OmnidirectionalAntennaElement,
                                                           elementRadius, numCircularElements,
                                                           elementBowRadius, elementBowAngle, numBowElements)

# plot antenna positions
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for antennaElement in antennaArray:
    xPos = antennaElement.positionVector[0]
    yPos = antennaElement.positionVector[1]
    zPos = antennaElement.positionVector[2]
    ax.scatter(xPos, yPos, zPos, marker='o', color='b')
plt.show()

# calculate weight vector
w = np.ones((len(antennaArray),1))

# calculate radiation pattern
poleResolution = 2**8
azimuthResolution = 2**8
poleAngles = np.array(range(poleResolution))/poleResolution*np.pi
azimuthAngles = np.array(range(azimuthResolution))/azimuthResolution*np.pi

arrayPattern = np.zeros((poleResolution, azimuthResolution), dtype=np.complex_)

for antennaIdx in range(len(antennaArray)):
    for thetaIdx, phiIdx in [(thetaIdx, phiIdx) for thetaIdx in range(poleResolution) 
                                                for phiIdx in range(azimuthResolution)]:
        theta = poleAngles[thetaIdx]
        phi = azimuthAngles[phiIdx]
    
        # assumption of a wave front (distance of rx to source -> infinite)
        # unit vector to the direction of wavefront
        directionVector = np.array([[np.sin(theta)*np.cos(phi)],
                                    [np.sin(theta)*np.sin(phi)],
                                    [np.cos(theta)]])
        # projection of element position vector to unit direction vector gives
        # the reference distance L
        antennaElement = antennaArray[antennaIdx]
        L = antennaElement.positionVector.reshape(1,-1) @ directionVector
        factor = np.squeeze(w[antennaIdx] * np.exp(1j*L*wavenumber))
        arrayPattern[thetaIdx,phiIdx] = arrayPattern[thetaIdx,phiIdx] + factor
        

arrayPatternLog = 10*np.log10(np.abs(np.real(arrayPattern)**2))
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
