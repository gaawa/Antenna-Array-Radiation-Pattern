from AntennaElement import AntennaElement

class OmnidirectionalAntennaElement(AntennaElement):
    def __init__(self, positionVector, orientationAng):
        super().__init__(positionVector, orientationAng)

    def getElementFactorElementBasis(signalAntElementBasis):
        return 1