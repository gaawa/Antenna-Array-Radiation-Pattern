from antenna_element.antenna_element_abs import AntennaElement
import numpy as np

class OmnidirectionalAntennaElement(AntennaElement):
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        super().__init__(positionVector, thetaOrientationVector, phiOrientationVector)

    def get_element_factor_element_basis(self, signalAngElementBasis):
        return np.ones(signalAngElementBasis.shape[1])