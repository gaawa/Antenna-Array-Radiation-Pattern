from antenna_element.antenna_element_abs import AntennaElement
import numpy as np

class HalfdirectionalAntennaElement(AntennaElement):
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        super().__init__(positionVector, thetaOrientationVector, phiOrientationVector)

   
    def get_element_factor_element_basis(self, signalAngElementBasis):
        # all angles with pole angle <= pi/2 gives 1, otherwise 0
        return 1*((signalAngElementBasis[0,:]).squeeze() <= np.pi/2)