from antenna_pattern.antenna_pattern_abs import AntennaPattern
import numpy as np

class HalfdirectionalAntennaPattern(AntennaPattern):
    """
    Defines an antenna element with hemispheric radiation pattern.
    """
    def __init__(self):
        super().__init__()

   
    def get_element_factor(self, signalVecLocal):
        # all angles with pole angle <= pi/2 gives 1, otherwise 0
        return 1*((signalVecLocal[0,:]).squeeze() <= np.pi/2)