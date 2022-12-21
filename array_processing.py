import matplotlib.pyplot as plt
from matplotlib import cm
from enumerations import PlotMode, BeamFormer
import numpy as np
from plotter import plot_radiation_pattern
from plotter import plot_antenna_topology
from utils import spherical_to_cartesian

def plot2d_generic(z_mat, x_vec, y_vec):
    """
    Generic 2D surf plot with x and y positions defined by the vector
    x_vec and y_vec

    Parameters
    ----------
    z_mat : 2D matrix
        Matrix with z values.
    x_vec : 1D vector
        Vector with x positions corresponding to values in z_mat.
    y_vec : 1D vector
        Vector with y positions corresponding to values in z_mat.

    Returns
    -------
    None.

    """
    iMesh, jMesh = np.meshgrid(x_vec, y_vec, indexing='ij')
    iMesh = iMesh/(2*np.pi)*360
    jMesh = jMesh/(2*np.pi)*360
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(jMesh, iMesh, z_mat, cmap=cm.jet)
    ax.set_xlabel('phi')
    ax.set_ylabel('theta')
    plt.show()
    

def calculate_radiation_pattern(antennaArray, w, poleAngles, azimuthAngles, wavenumber):
    """
    Calculates the array radiation pattern by calculating the signal amplification
    of the signal coming from the directions defined by poleAngles and azimuthAngles.

    Parameters
    ----------
    antennaArray : AntennaArray
        AntennaArray object that defines the antenna array.
    w : 1D vector
        Weighting factor for the array element in AntennaArray.arrayElements.
    poleAngles : 1D vector
        The pole angles to be scanned.
    azimuthAngles : 1D vector
        The azimuth angles to be scanned.
    wavenumber : floating point
        2*pi/lambda of the system.

    Returns
    -------
    arrayPattern : 2D matrix
        Matrix where its values represent received signal amplitude where the
        axis0 and axis1 of the matrix corresponds to the direction of arrival 
        given by the poleAngles and azimuthAngles respectively.

    """
    poleResolution = poleAngles.size
    azimuthResolution = azimuthAngles.size

    arrayPattern = np.zeros((poleResolution, azimuthResolution), dtype=np.complex_)

    # create cartesian direction vectors as a function of 
    # spherical coordinate parameters
    spherical_param_mat = np.array([(theta, phi) for theta in poleAngles 
                                                 for phi in azimuthAngles]).T
    thetaVec = spherical_param_mat[0,:][:]
    phiVec = spherical_param_mat[1,:][:]
    dv_matrix = spherical_to_cartesian(1, thetaVec, phiVec)

    # create list of antenna position vectors from all elements
    pv_list = [ antennaElement.positionVector for antennaElement in antennaArray ]
    pv_matrix = np.array(pv_list).T

    # projection of direction to position vector gives the propagation distance
    # relative to reference point (0,0,0).
    L_matrix = pv_matrix.T @ dv_matrix

    # get antenna element factor
    G_ant = [ antennaElement.get_element_factor_array_basis(dv_matrix)
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
    """
    Performs the analysis of the antenna array.
    1. Plots position and orientation of each antenna elements of the array
    2. Calculates the weighting factors with the method given by 'beamformer'
    3. Calculates the radiation pattern
    4. Plots the radiation pattern
    5. Print out amplitude and phase shift of the weighting factors

    Parameters
    ----------
    antennaArray : AntennaArray
        Instance of AntennaArray object that defines the array.
    wavenumber : floating point
        2*pi/lambda of the system.
    beamformer : enumerations.BeamFormer or None
        The type of beamforming to perform.
    beamformerTheta : floating point
        Pole angle of the beamforming direction.
    beamformerPhi : floating point
        Azimuth angle of the beamforming direction.
    poleResolution : integer
        How many angle steps the pole angle has for plotting.
    azimuthResolution : integer
        How man angle steps the azimuth angle has for plotting.
    plotThetaRange : floating poing
        The angle range of the pole angle for plotting.
    plotPhiRange : floating poing
        The angle range of the azimuth angle for plotting.
    plotMode : enumerations.PlotMode
        Plotting mode for radiation pattern plot.

    Returns
    -------
    None.

    """
    arrayElements = antennaArray.arrayElements
    ## plot antenna positions
    # plot_antenna_topology.plotDots(antennaArray)
    plot_antenna_topology.plot_quiver(arrayElements, 
                                     ax_size=antennaArray.ax_size, 
                                     arrow_size=antennaArray.arrow_size)
    
    ## calculate weight vector (Beamforming)
    if not beamformer:
        w = np.ones((len(arrayElements),1))
    elif beamformer is BeamFormer.Projection:
        # create cartesian direction vector in the direction of desired beam steering
        dv = spherical_to_cartesian(1, 
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
    arrayPattern = calculate_radiation_pattern(arrayElements, w, poleAngles, azimuthAngles,
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
        
