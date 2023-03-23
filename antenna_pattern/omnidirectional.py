from antenna_pattern.antenna_pattern_abs import AntennaPattern
import numpy as np

class OmnidirectionalAntennaPattern(AntennaPattern):
    """
    Defines omnidirectional antenna array.
    """
    def __init__(self):
        super().__init__()

    def get_element_factor(self, signalVecLocal):
        return np.ones(signalVecLocal.shape[1])