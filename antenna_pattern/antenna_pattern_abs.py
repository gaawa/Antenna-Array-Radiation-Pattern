from abc import ABC, abstractmethod
import numpy as np

class AntennaPattern(ABC):
    """
    Abstract base class for all antenna pattern classes

    All child classes must implement get_element_factor function,
    which defines the antenna element radiation pattern
    """
    
    def __init__(self):
        super().__init__()

    @abstractmethod
    def get_element_factor(self, signalVecLocal):
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
            'self.phiOrientationVector' (in AntennaElement object).

        Returns
        -------
        1D vector
            signal amplification due to the antenna radiation pattern corressponding
            to the direction of arrival given by 'signalVecLocal'.

        """
        pass