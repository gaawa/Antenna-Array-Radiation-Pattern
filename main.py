import matplotlib.pyplot as plt
from matplotlib import cm
from ArrayGenerator import LinearArray, CircularArray, LinearCircularArray, BowCircularArray
from enumerations import ArrayFormation, PlotMode, BeamFormer
import numpy as np
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
    

def calculateRadiationPattern(antennaArray, w, poleAngles, azimuthAngles,
                              wavenumber):
    poleResolution = poleAngles.size
    azimuthResolution = azimuthAngles.size

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
    # plot2d_generic(G_ant[0,:].reshape((poleResolution, azimuthResolution)), poleAngles, azimuthAngles)

    # phase shift of the signal is then the relative distance multiplied by the wave number
    # calculate signal amplitude A(i,j)
    # axis 0: at the antenna element i
    # axis 1: of signal direction j
    factorMatrix = w * G_ant * np.exp(1j*L_matrix*wavenumber) # todo: antenna element factor (map to function)

    arrayFactorVector = np.sum(factorMatrix, axis=0)
    arrayPattern = arrayFactorVector.reshape((poleResolution, azimuthResolution))
    
    return arrayPattern


def array_processing(antennaArray, wavenumber, 
                     beamformer, beamformerTheta, beamformerPhi,
                     poleResolution, azimuthResolution,
                     plotThetaRange, plotPhiRange,
                     plotMode):
    arrayElements = antennaArray.arrayElements
    ## plot antenna positions
    # plot_antenna_topology.plotDots(antennaArray)
    plot_antenna_topology.plotQuiver(arrayElements, 
                                     ax_size=antennaArray.ax_size, 
                                     arrow_size=antennaArray.arrow_size)
    
    ## calculate weight vector (Beamforming)
    if not beamformer:
        w = np.ones((len(arrayElements),1))
    elif beamformer is BeamFormer.Projection:
        # create cartesian direction vector in the direction of desired beam steering
        dv = coordinate_utils.spherical_to_cartesian(1, 
                                                     beamformerTheta, 
                                                     beamformerPhi).reshape(-1,1)
    
        # create list of antenna position vectors from all elements
        pv_list = [ antennaElement.positionVector for antennaElement in arrayElements ]
        pv_matrix = np.array(pv_list).T
        
        # projection of direction to position vector gives the propagation distance
        # relative to reference point (0,0,0).
        L_matrix = pv_matrix.T @ dv
        
        deltaPhase = L_matrix*wavenumber
        w = np.exp(-1j*deltaPhase)
    
    ## calculate radiation pattern
    poleAngles = np.array(range(poleResolution))/(poleResolution-1)*plotThetaRange
    azimuthAngles = np.array(range(azimuthResolution))/(azimuthResolution-1)*plotPhiRange
    arrayPattern = calculateRadiationPattern(arrayElements, w, poleAngles, azimuthAngles,
                                             wavenumber)
    
    if plotMode is PlotMode.plot2d:    
        plot_radiation_pattern.plot2d(arrayPattern, poleAngles, azimuthAngles)
    else:
        dv = plot_radiation_pattern.plot3d(arrayPattern, 
                                           poleAngles, 
                                           azimuthAngles)
        
    ## print out weighting factors
    w_amp = np.abs(w)
    w_deg = (np.angle(w))*360/(2*np.pi)
    
    for i in range(w_amp.size):
        print("Antenna {}: {} amp {} deg".format(i, w_amp[i], w_deg[i]))
        
