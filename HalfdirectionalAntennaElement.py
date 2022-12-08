from AntennaElement import AntennaElement
import numpy as np

class HalfdirectionalAntennaElement(AntennaElement):
    def __init__(self, positionVector, orientationAng):
        super().__init__(positionVector, orientationAng)

    def getElementFactorElementBasis(signalAntElementBasis):
        if signalAntElementBasis[0] <= np.pi:
            return 1
        return 0