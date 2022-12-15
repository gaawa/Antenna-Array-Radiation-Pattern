from AntennaElement import AntennaElement
import numpy as np

class HalfdirectionalAntennaElement(AntennaElement):
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        super().__init__(positionVector, thetaOrientationVector, phiOrientationVector)

    def getElementFactorElementBasis(self, signalAngElementBasis):
        # all angles with pole angle <= pi/2 gives 1, otherwise 0
        return 1*((signalAngElementBasis[0,:]).squeeze() <= np.pi/2)