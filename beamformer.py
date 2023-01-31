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

def synthesis_beamformer(antennaArray, beamformerPropTuples, wavenumber):
    """
    Creates a beamformer according to properties defined by 'beamformerPropTuples'.
    The properties specifies the array factor of the beamformer at a given 
    signal direction.
    The beamformer is calculated by solving a linear equation system.

    Parameters
    ----------
    antennaArray : AntennaArray
        AntennaArray object.
    beamformerPropTuples : List[Tuple(antennaFactor, theta, phi)]
        List of tuples that contains the following beamformer properties.
        
        antennaFactor : floating point
            Signal amplification for signals from the direction defined by
            the angle theta and phi.
        theta : floating point
            Pole angle of the directiokn of arrival of the signal subject to
            the amplification defined by 'antennaFactor'
        phi : floating point
            Azimuth angle of the directiokn of arrival of the signal subject to
            the amplification defined by 'antennaFactor'
    wavenumber : floating point
        2*pi/lambda of the system.

    Returns
    -------
    w : vector in 2D matrix (Nx1 matrix)
        Weighting factor corresponding to antenna elements in 
        antennaArray.antennaElements

    """
    antennaFactors = np.array([propTuple[0] for propTuple in beamformerPropTuples]).reshape(-1,1)
    thetas = np.array([propTuple[1] for propTuple in beamformerPropTuples])
    phis = np.array([propTuple[2] for propTuple in beamformerPropTuples])
    
    # create cartesian direction vector in the direction of desired beam steering
    directionVecs = spherical_to_cartesian(1, thetas, phis)
    
    # create list of antenna position vectors from all elements
    positionVecs = np.array([ antennaElement.positionVector 
                             for antennaElement in antennaArray.arrayElements ]).T
    
    # projection of direction to position vectors gives the propagation distances
    # relative to reference point (0,0,0).
    # the index is
    # axis 0: of signal direction i
    # axis 1: antenna element j
    L_matrix = directionVecs.T @ positionVecs
    
    # get antenna element factor (antenna gain)
    # axis 0: of signal direction i
    # axis 1: antenna element j
    G_ant = [ antennaElement.get_element_factor_array_basis(directionVecs)
                 for antennaElement in antennaArray.arrayElements]
    G_ant = np.array(G_ant).T
    
    # A-matrix where the index is
    # axis 0: of signal direction i
    # axis 1: at the antenna element j
    A = G_ant*np.exp(1j*L_matrix*wavenumber)
    # w, res, rnk, s = np.linalg.lstsq(A, antennaFactors)
    
    # get pseudo inverse of A
    A_inv = A.conjugate().T @ np.linalg.inv(A@A.conjugate().T)
    
    w = A_inv @ antennaFactors
    return w
    
def synthesis_canceller_beamformer(antennaArray, w_q, beamformerPropTuples, wavenumber):
    """
    Creates a beamformer according to properties defined by 'beamformerPropTuples'.
    The properties specifies the array factor of the beamformer at a given 
    signal direction.
    The beamformer is calculated by combining a quiescent beam with a second
    'cancelling' beam that substracts signal power from specified directions.

    Parameters
    ----------
    antennaArray : AntennaArray
        AntennaArray object.
    beamformerPropTuples : List[Tuple(antennaFactor, theta, phi)]
        List of tuples that contains the following beamformer properties.
        
        antennaFactor : floating point
            Signal amplification for signals from the direction defined by
            the angle theta and phi.
        w_q : 1D/2D vector
            Beamforming weights for the quiescent beam.
        theta : floating point
            Pole angle of the directiokn of arrival of the signal subject to
            the amplification defined by 'antennaFactor'
        phi : floating point
            Azimuth angle of the directiokn of arrival of the signal subject to
            the amplification defined by 'antennaFactor'
    wavenumber : floating point
        2*pi/lambda of the system.

    Returns
    -------
    w : vector in 2D matrix (Nx1 matrix)
        Weighting factor corresponding to antenna elements in 
        antennaArray.antennaElements

    """
    # desired antenna factors and its corresponding directions
    AF_d = np.array([propTuple[0] for propTuple in beamformerPropTuples]).reshape(-1,1)
    thetas = np.array([propTuple[1] for propTuple in beamformerPropTuples])
    phis = np.array([propTuple[2] for propTuple in beamformerPropTuples])
    w_q = w_q.reshape(-1,1)
    
    # create cartesian direction vector in the direction of specified
    # desired beam steering array factor
    directionVecs = spherical_to_cartesian(1, thetas, phis)
    
    # create list of antenna position vectors from all elements
    positionVecs = np.array([ antennaElement.positionVector 
                             for antennaElement in antennaArray.arrayElements ]).T
    
    # projection of direction to position vectors gives the propagation distances
    # (the delta phase) relative to reference point (0,0,0).
    # deltaPhases matrix where the index is
    # axis 0: of signal direction i
    # axis 1: antenna element j
    L_matrix = directionVecs.T @ positionVecs
    deltaPhases = L_matrix*wavenumber
    
    # get antenna element factor
    G_ant = [ antennaElement.get_element_factor_array_basis(directionVecs)
                 for antennaElement in antennaArray.arrayElements]
    G_ant = np.array(G_ant).T
    
    # A-matrix where the index is
    # axis 0: at the signal direction i
    # axis 1: at the antenna element j
    A = G_ant*np.exp(1j*deltaPhases)
    
    # calculate quiescent antenna factor
    AF_q = A @ w_q
    
    # get pseudo inverse of A-matrix
    A_inv = A.conjugate().T @ np.linalg.inv(A@A.conjugate().T)
    
    # calculate cancellation beam weights
    w_c = A_inv @ (AF_d - AF_q)
    
    return w_c + w_q