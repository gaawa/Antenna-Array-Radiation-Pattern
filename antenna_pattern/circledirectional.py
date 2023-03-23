from antenna_pattern.antenna_pattern_abs import AntennaPattern
import numpy as np

class CircledirectionalAntennaPattern(AntennaPattern):
    """
    Defines an antenna element with circle directional radiation pattern.
    This means that the radiation pattern has a form of a displaced circle.
    The circle is however in db scale.
    The scale range is from dBmin to dBmax.
    The default is dBmin=-20 and dBmax=0.
    This means that on the front side, the gain is 0dB, on the back side -20dB.
    Values are converted to linear scale for output
    """
    def __init__(self, dBmin=-20, dBmax=0):
        super().__init__()
        self.dBmin = dBmin 
        self.dBmax = dBmax
   
    def get_element_factor(self, signalVecLocal):
        # all angles with pole angle <= pi/2 gives 1, otherwise 0
        thetas = signalVecLocal[0,:].squeeze()
        rdB = (self.dBmax-self.dBmin)*np.cos(thetas) + self.dBmin
        rlin = 10**(rdB/20)
        return rlin