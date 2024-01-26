from array_processing import array_processing
from antenna_array import LinearArray, CircularArray, LinearCircularArray, BowCircularArray
from enumerations import PlotMode, BeamFormer
import scipy.constants as spc
import numpy as np
from antenna_pattern.omnidirectional import OmnidirectionalAntennaPattern
from antenna_pattern.halfdirectional import HalfdirectionalAntennaPattern
from antenna_pattern.circledirectional import CircledirectionalAntennaPattern
from antenna_pattern.csv_file_pattern import CsvFilePattern
from beamformer import no_beamformer, projection_beamformer, partial_projection_beamformer
from beamformer import synthesis_beamformer, synthesis_canceller_beamformer
import matplotlib.pyplot as plt

## set physical constants
freq = 3e9
wavelength = spc.speed_of_light/freq
wavenumber = 2*np.pi/wavelength

antennaPattern = CsvFilePattern('antenna_pattern/Gain_lin_simulation_johann.csv', debug=False, fastMode=False)
antennaPattern = CsvFilePattern('antenna_pattern/Thermocomp_0-45.csv', debug=False, fastMode=True)
# antennaPattern = CircledirectionalAntennaPattern()

angResolution = 180
azimuthOffset = np.pi/4
plotThetaRange = np.pi

poleAngles = np.array(range(angResolution))/(angResolution-1)*plotThetaRange

azimuthAnglesZero = np.zeros(angResolution)+azimuthOffset
azimuthAnglesPi = np.pi*np.ones(angResolution)+azimuthOffset



angles = np.hstack((np.vstack((np.flip(poleAngles), azimuthAnglesZero)), 
                         np.vstack((poleAngles, azimuthAnglesPi))))
elementFactor = antennaPattern.get_element_factor(angles)

fig = plt.figure(figsize=[10, 10])
ax = fig.add_subplot()
ax.plot(10*np.log10(elementFactor), 'x')
plt.show()