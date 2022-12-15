from abc import ABC, abstractmethod
import numpy as np

class AntennaElement(ABC):
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        """
        

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
    def getElementFactorElementBasis(self, signalVecLocal):
        # the implementation must imply that the input is 2D column vecotors.
        # it must cover the case when the input is 2xN matrix,
        # where the axis 0 is [theta, phi] and axis 1 the multiple different
        # angle directions.
        # The output must then be an array, each element corresponding to the
        # direction in the axis 1
        pass

    def getElementFactorArrayBasis(self, signalVecGlobal):
        # convert input to 3x1 column vector or 3xN matrix of column vectors
        signalVecGlobal = signalVecGlobal.reshape(3,-1)
        
        # calculate local poleAngle (theta) and Azimuth (phi)
        theta = np.arccos( np.dot(self.thetaOrientationVector, signalVecGlobal) /
                           ( np.linalg.norm(self.thetaOrientationVector)
                            *np.linalg.norm(signalVecGlobal, axis=0) ))
        
        signalThetaAxisProjection = (np.dot(self.thetaOrientationVector, signalVecGlobal) 
                                    / np.linalg.norm(self.thetaOrientationVector)**2
                                    * self.thetaOrientationVector.reshape((3,-1)) )                       
        signalAntPlaneProjection = signalVecGlobal - signalThetaAxisProjection
        phi = np.arccos( np.dot(self.phiOrientationVector, signalAntPlaneProjection) /
                           ( np.linalg.norm(self.phiOrientationVector)
                            *np.linalg.norm(signalAntPlaneProjection, axis=0) ) )
        
        signalAngLocal = np.array([ theta,
                                    phi ])
        
        return self.getElementFactorElementBasis(signalAngLocal)