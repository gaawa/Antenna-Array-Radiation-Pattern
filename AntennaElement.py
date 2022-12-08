from abc import ABC, abstractmethod

class AntennaElement(ABC):
    def __init__(self, positionVector, orientationAng):
        self.positionVector = positionVector
        self.orientationAng = orientationAng
        super().__init__()

    @abstractmethod
    def getElementFactorElementBasis(signalAngElementBasis):
        pass

    def getElementFactorArrayBasis(self, signalAngArrayBasis):
        signalAngElementBasis = signalAngArrayBasis - self.orientationAng
        return self.getElementFactorArrayBasis(signalAngElementBasis)