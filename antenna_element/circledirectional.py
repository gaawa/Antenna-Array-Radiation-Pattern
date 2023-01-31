from antenna_element.antenna_element_abs import AntennaElement
import numpy as np

class CircledirectionalAntennaElement(AntennaElement):
    """
    Defines an antenna element with circle directional radiation pattern.
    This means that the radiation pattern has a form of a displaced circle.
    The circle is however in db scale.
    The scale range is from -20 to 0 dB.
    On the front side, the gain is 0dB, on the back side -10dB.
    Values are converted to linear scale for output
    """
    def __init__(self, positionVector, thetaOrientationVector, phiOrientationVector):
        super().__init__(positionVector, thetaOrientationVector, phiOrientationVector)
        self.dBmin = -20
        self.dBmax = 0
   
    def get_element_factor_element_basis(self, signalAngElementBasis):
        # all angles with pole angle <= pi/2 gives 1, otherwise 0
        thetas = signalAngElementBasis[0,:].squeeze()
        rdB = (self.dBmax-self.dBmin)*np.cos(thetas) + self.dBmin
        rlin = 10**(rdB/20)
        return rlin