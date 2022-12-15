from AntennaElement import AntennaElement
import numpy as np

class OmnidirectionalAntennaElement(AntennaElement):
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        super().__init__(positionVector, thetaOrientationVector, phiOrientationVector)

    def getElementFactorElementBasis(self, signalAngElementBasis):
        return np.ones(signalAngElementBasis.shape[1])