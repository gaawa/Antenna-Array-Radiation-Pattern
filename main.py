import matplotlib.pyplot as plt
from matplotlib import cm
import ArrayGenerator
from OmnidirectionalAntennaElement import OmnidirectionalAntennaElement
from HalfdirectionalAntennaElement import HalfdirectionalAntennaElement
import scipy.constants as spc
import numpy as np
from enum import Enum, auto
import equidist_pts_sphere
import coordinate_utils
import plot_radiation_pattern
import plot_antenna_topology

def plot2d_generic(R, poleAngles, azimuthAngles):
    iMesh, jMesh = np.meshgrid(poleAngles, azimuthAngles, indexing='ij')
    iMesh = iMesh/(2*np.pi)*360
    jMesh = jMesh/(2*np.pi)*360
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(jMesh, iMesh, R, cmap=cm.jet)
    ax.set_xlabel('phi')
    ax.set_ylabel('theta')
    plt.show()
    
class ArrayFormation(Enum):
    Linear = auto()
    Circular = auto()
    LinearCircular = auto()
    BowCircular = auto()
    
class PlotMode(Enum):
    plot2d = auto()
    plot3d = auto()

# spawn antenna elements in linear array
freq = 3e9
wavelength = spc.speed_of_light/freq
wavenumber = 2*np.pi/wavelength
elementDistance = wavelength/2
arrayFormation = ArrayFormation.BowCircular
AntennaElementClass = HalfdirectionalAntennaElement
plotMode = PlotMode.plot3d

plotThetaRange = np.pi
plotPhiRange = 2*np.pi
poleResolution = 2**6
azimuthResolution = 2**6

if arrayFormation is ArrayFormation.Linear:
    numElements = 10
    ax_size = 1.2*numElements*elementDistance
    arrow_size = 0.5*elementDistance
    antennaArray = ArrayGenerator.generateLinearArray(AntennaElementClass,
                                                      elementDistance, numElements)


elif arrayFormation is ArrayFormation.Circular:
    numElements = 10
    elementRadius = wavelength/(4*np.sin(np.pi/numElements))
    poleAng = np.pi/2
    ax_size = 1.2*2*elementRadius
    arrow_size = 0.5*(2*np.pi*elementRadius/numElements)
    antennaArray = ArrayGenerator.generateCircularArray(AntennaElementClass,
                                                        elementRadius, numElements,
                                                        poleAng=poleAng)


elif arrayFormation is ArrayFormation.LinearCircular:
    numCircularElements = 10
    elementRadius = wavelength/(4*np.sin(np.pi/numCircularElements))
    numLinearElements = 4
    elementLinearDistance = wavelength/2
    ax_size = 1.2 * max(numLinearElements*elementLinearDistance,
                    2*elementRadius)
    arrow_size = 0.5 * min(elementDistance,
                           2*np.pi*elementRadius/numCircularElements)
    antennaArray = ArrayGenerator.generateLinearCircularArray(AntennaElementClass,
                                                              elementRadius, numCircularElements,
                                                              elementLinearDistance, numLinearElements)


elif arrayFormation is ArrayFormation.BowCircular:
    numCircularElements = 10
    elementRadius = wavelength/(4*np.sin(np.pi/numCircularElements))
    numBowElements = 4
    elementBowAngle = np.pi/2
    ax_size = 1.2 * 2*elementRadius
    arrow_size = 0.5 * 2*np.pi*elementRadius/numCircularElements
    elementBowRadius = elementDistance/(2*np.sin(elementBowAngle/(2*numBowElements)))
    antennaArray = ArrayGenerator.generateBowCircularArray(AntennaElementClass,
                                                           elementRadius, numCircularElements,
                                                           elementBowRadius, elementBowAngle, numBowElements)


# plot antenna positions
# plot_antenna_topology.plotDots(antennaArray)
plot_antenna_topology.plotQuiver(antennaArray, ax_size=ax_size, arrow_size=arrow_size)

# calculate weight vector (Beamforming)
w = np.ones((len(antennaArray),1))

# calculate radiation pattern
poleAngles = np.array(range(poleResolution))/(poleResolution-1)*plotThetaRange
azimuthAngles = np.array(range(azimuthResolution))/(azimuthResolution-1)*plotPhiRange

arrayPattern = np.zeros((poleResolution, azimuthResolution), dtype=np.complex_)

# create cartesian direction vectors as a function of 
# spherical coordinate parameters
spherical_param_mat = np.array([(theta, phi) for theta in poleAngles 
                                             for phi in azimuthAngles]).T
thetaVec = spherical_param_mat[0,:][:]
phiVec = spherical_param_mat[1,:][:]
dv_matrix = coordinate_utils.spherical_to_cartesian(1, thetaVec, phiVec)

# create list of antenna position vectors from all elements
pv_list = [ antennaElement.positionVector for antennaElement in antennaArray ]
pv_matrix = np.array(pv_list).T

# projection of direction to position vector gives the propagation distance
# relative to reference point (0,0,0).
L_matrix = pv_matrix.T @ dv_matrix

# get antenna element factor
G_ant = [ antennaElement.getElementFactorArrayBasis(dv_matrix)
             for antennaElement in antennaArray]
G_ant = np.array(G_ant)
plot2d_generic(G_ant[0,:].reshape((poleResolution, azimuthResolution)), poleAngles, azimuthAngles)

# phase shift of the signal is then the relative distance multiplied by the wave number
# calculate signal amplitude A(i,j)
# axis 0: at the antenna element i
# axis 1: of signal direction j
factorMatrix = w * G_ant * np.exp(1j*L_matrix*wavenumber) # todo: antenna element factor (map to function)

arrayFactorVector = np.sum(factorMatrix, axis=0)
arrayPattern = arrayFactorVector.reshape((poleResolution, azimuthResolution))

if plotMode is PlotMode.plot2d:    
    plot_radiation_pattern.plot2d(arrayPattern, poleAngles, azimuthAngles)
else:
    plot_radiation_pattern.plot3d(arrayPattern, poleAngles, azimuthAngles)
