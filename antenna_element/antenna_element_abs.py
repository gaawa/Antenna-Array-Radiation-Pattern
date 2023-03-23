from abc import ABC, abstractmethod
import numpy as np

class AntennaElement(ABC):
    """
    Abstract base class for all antenna element classes.
    
    All child classes must implement get_element_factor_element_basis function,
    which defines the antenna element radiation pattern.
    """
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        """
        This ABC requires the coordinate position of the antenna element 'positionVector'
        and two vectors 'thetaOrientationVector' and 'phiOrientationVector' which
        defines the axis of 0 pole angle and 0 azimuth angle respectively.

        Parameters
        ----------
        positionVector : 3 numpyarray
            position of antenna element relative to the global origin.
        thetaOrientationVector : 3 numpyarray
            global coordinate unit vector pointing to the direction where 
            the local antenna theta axis points.
        phiOrientationVector : 3 numpyarray
            global coordinate unit vector pointing to the direction where 
            the local antenna phi axis points.

        Returns
        -------
        None.

        """
        # vectors in the format np.array([x,y,z])
        self.positionVector = positionVector.reshape(3)
        self.thetaOrientationVector = thetaOrientationVector.reshape(3)
        self.phiOrientationVector = phiOrientationVector.reshape(3)
        super().__init__()

    @abstractmethod
    def get_element_factor_element_basis(self, signalVecLocal):
        """
        Defines the radiation pattern of the antenna element.
        Implementation returns signal amplification of the signal with 
        direction of arrival defined by signalVecLocal (local coordinate system), 
        due to the antenna radiation pattern.
        
        The implementation must imply that 'signalVecLocal' is 2D column vecotors.
        It must cover the case where the input is 2xN matrix
        whereby the axis 0 is [theta, phi] and axis 1 the multiple different
        angle directions.
        The output must then be an array, each element corresponding to the
        direction in the axis 1

        Parameters
        ----------
        signalVecLocal : 1D or 2D matrix
            vector or matrix containing one or multiple vectors of [theta, phi],
            which is the pole angle and azimuth angle respectively.
            For multiple vectors, the dimension shall be 2xN.
            These vectors are unit vectors in spherical coordinate system where
            the basis vectors are defined by 'self.thetaOrientationVector' and
            'self.phiOrientationVector'.

        Returns
        -------
        1D vector
            signal amplification due to the antenna radiation pattern corressponding
            to the direction of arrival given by 'signalVecLocal'.

        """
        pass
    
    def get_element_factor_array_basis(self, signalVecGlobal):
        """
        Calculates the signal amplification due to the antenna radiation pattern
        for direction of arrival defined by 'signalVecGlobal' (global coordinate system).
        
        Converts the input vectors, which are in cartesian coordinate system with
        global basis vectors to unit vectors in spherical coordinate system with 
        local basis vectors defined by 'self.thetaOrientationVector' and
        'self.phiOrientationVector'.
        
        Parameters
        ----------
        signalVecGlobal : 1D vector or 2D matrix
            vector or matrix containing one or multiple vectors of [x,y,z] 
            direction of arrival vector.
            For multiple vectors, the dimension shall be 3xN.
            

        Returns
        -------
        1D vector
            signal amplification due to the antenna radiation pattern 
            corresponding to the direction of arrival given by 'signalVecGlobal'.

        """
        # convert input to 3x1 column vector or 3xN matrix of column vectors
        signalVecGlobal = signalVecGlobal.reshape(3,-1)
        
        # calculate local poleAngle (theta)
        # t.u' = |t|*|u|*cos(theta) -> solve to theta
        theta = np.arccos( np.dot(self.thetaOrientationVector, signalVecGlobal) /
                           (  np.linalg.norm(self.thetaOrientationVector)
                            * np.linalg.norm(signalVecGlobal, axis=0) ))
        
        # calculate local Azimuth (phi)
        # Step 1) Project to theta-axis: u_t = proj u->t = u.t' / (||t||^2) * t
        signalThetaAxisProjection = (np.dot(self.thetaOrientationVector, signalVecGlobal) 
                                    / np.linalg.norm(self.thetaOrientationVector)**2
                                    * self.thetaOrientationVector.reshape((3,-1)) )
        # Step 2) Obtain signal vector projected onto antenna plane: u_Ap = proj u->Ap = u - u_t
        signalAntPlaneProjection = signalVecGlobal - signalThetaAxisProjection
        # Step 3) Obtain cos(phi) between u_Ap and phi-axis: p.u_Ap' = ||p|| * ||u_Ap|| * cos(phi)
        cosphi = (np.dot(self.phiOrientationVector, signalAntPlaneProjection) /
                    (  np.linalg.norm(self.phiOrientationVector)
                     * np.linalg.norm(signalAntPlaneProjection, axis=0) ))
        # Step 4) Obtain sin(phi) between u_Ap and phi-axis: p x u_Ap = ||p|| * ||u_A|| * sin(phi) * n
        # Get value along normal vector n by projecting (p x u_Ap)/(||p|| * ||u_A||) onto n.
        # Note that normal vector n is equivalent to theta-axis unit vector
        sinphi_normal = ( np.cross(self.phiOrientationVector[:,None], signalAntPlaneProjection, axis=0) /
                            (  np.linalg.norm(self.phiOrientationVector)
                            * np.linalg.norm(signalAntPlaneProjection, axis=0) ))
        sinphi = np.dot(self.thetaOrientationVector, sinphi_normal)
        # Step 5) Calculate full [0, 2*pi] phi angle
        phi = np.arctan2(sinphi, cosphi)
        phi = np.nan_to_num(phi)

        signalAngLocal = np.array([ theta,
                                    np.nan_to_num(phi, nan=0.0) ])
        return self.get_element_factor_element_basis(signalAngLocal)