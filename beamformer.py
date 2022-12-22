from abc import ABC, abstractmethod
from utils import spherical_to_cartesian
import numpy as np

def no_beamformer(antennaArray):
    """
    No beamforming. All weighting factors are set to 1.

    Parameters
    ----------
    antennaArray : AntennaArray
        AntennaArray object.

    Returns
    -------
    w : vector in 2D matrix (Nx1 matrix)
        Weighting factor corresponding to antenna elements in 
        antennaArray.antennaElements

    """
    numElements = len(antennaArray.arrayElements)
    w = np.ones((numElements, 1))
    return w

def projection_beamformer(antennaArray, beamformerTheta, beamformerPhi, wavenumber):
    """
    Simple beamformer obtained by phase shift calculations assuming a reception
    of a wave front from the direction defined by beamformerTheta and beamformerPhi.
    Phase shift can be calculated with the wavenumber and the projection of the
    direction of arrival vector onto the antenna position vector.

    Parameters
    ----------
    antennaArray : AntennaArray
        AntennaArray object.
    beamformerTheta : floating point
        Pole angle of the direction of arrival of the desired signal.
    beamformerPhi : floating point
        Azimuth angle of the direction of arrival of the desired signal.
    wavenumber : floating point
        2*pi/lambda of the system.

    Returns
    -------
    w : vector in 2D matrix (Nx1 matrix)
        Weighting factor corresponding to antenna elements in 
        antennaArray.antennaElements

    """    
    # create cartesian direction vector in the direction of desired beam steering
    dv = spherical_to_cartesian(1, 
                                beamformerTheta, 
                                beamformerPhi).reshape(-1,1)

    # create list of antenna position vectors from all elements
    pv_list = [ antennaElement.positionVector for antennaElement in antennaArray.arrayElements ]
    pv_matrix = np.array(pv_list).T
    
    # projection of direction to position vector gives the propagation distance
    # relative to reference point (0,0,0).
    L_matrix = pv_matrix.T @ dv
    
    deltaPhase = L_matrix*wavenumber
    w = np.exp(-1j*deltaPhase)
    
    return w

def partial_projection_beamformer(antennaArray, beamformerTheta, beamformerPhi, 
                                  wavenumber, active_angle=np.pi/2):
    """
    Same as the function projection_beamformer but only with antennas that
    are oriented within the angle of 'active_angle' relative to the direction
    of the desired signal. Weighting factors for all other antenna elements
    are set to 0.

    Parameters
    ----------
    antennaArray : AntennaArray
        AntennaArray object.
    beamformerTheta : floating point
        Pole angle of the direction of arrival of the desired signal.
    beamformerPhi : floating point
        Azimuth angle of the direction of arrival of the desired signal.
    wavenumber : floating point
        2*pi/lambda of the system.
    active_angle : floating point, optional
        The angle relative from the beamforming direction.
        Only antenna arrays oriented within this angle are activated. 
        The default is np.pi/2.

    Returns
    -------
    w : TYPE
        DESCRIPTION.

    """
    # create cartesian direction vector in the direction of desired beam steering
    dv = spherical_to_cartesian(1, 
                                beamformerTheta, 
                                beamformerPhi).reshape(-1,1)

    # create list of antenna position vectors from all elements
    pv_list = [ antennaElement.positionVector for antennaElement in antennaArray.arrayElements ]
    pv_matrix = np.array(pv_list).T
    
    # projection of direction to position vector gives the propagation distance
    # relative to reference point (0,0,0).
    L_matrix = pv_matrix.T @ dv
    
    deltaPhase = L_matrix*wavenumber
    w = np.exp(-1j*deltaPhase)
    
    # mask out weighting factors for antennas that are oriented >active_angle away
    # from the beamforming direction vector
    ov_list = [ antennaElement.thetaOrientationVector for antennaElement in antennaArray.arrayElements]
    ov_matrix = np.array(ov_list).T
    od_dot = ov_matrix.T @ dv
    angle_delta = np.arccos(od_dot/(np.linalg.norm(ov_matrix, axis=0).reshape(od_dot.shape)
                                   *np.linalg.norm(dv)))
    mask = angle_delta <= active_angle
    w = w * mask
    
    return w