# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 15:19:42 2022

@author: ThinkWin
"""
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

## set physical constants
freq = 3e9
wavelength = spc.speed_of_light/freq
wavenumber = 2*np.pi/wavelength

## select antenna pattern type
antennaPattern = HalfdirectionalAntennaPattern()
antennaPattern = OmnidirectionalAntennaPattern()
antennaPattern = CircledirectionalAntennaPattern()
antennaPattern = CsvFilePattern('antenna_pattern/Gain_lin_simulation.csv', debug=False, fastMode=True)

## spawn antenna elements with specified array topology
# Single antenna element for testing
antennaArray = LinearArray(antennaPattern, wavelength, numElements=1, elementDistanceFactor=0.5)
# Linear array with 8 elements
antennaArray = LinearArray(antennaPattern, wavelength, numElements=4, elementDistanceFactor=0.5)
# Spherical conformal array
antennaArray = BowCircularArray(antennaPattern, wavelength, circularAng=np.pi/4, bowAng=np.pi/2)
# Zylindrical conformal array
# antennaArray = LinearCircularArray(antennaPattern, wavelength) 

## set desired beamforming direction
beamformerTheta = np.pi/4
beamformerPhi = np.pi/4

## calculate weighting factors
w = no_beamformer(antennaArray)
# w = projection_beamformer(antennaArray, beamformerTheta, beamformerPhi, wavenumber)
# w = partial_projection_beamformer(antennaArray, beamformerTheta, beamformerPhi, 
#                                   wavenumber, active_angle=np.pi/2)
# w = synthesis_beamformer(antennaArray, [(6, beamformerTheta, beamformerPhi),
#                                         (0, np.radians(30), np.radians(250)),
#                                         (0, np.radians(90), np.radians(50)),
#                                         (0, np.radians(50), np.radians(145)),
#                                         (0, np.radians(50), np.radians(320))], 
#                                          wavenumber)
# w = synthesis_canceller_beamformer(antennaArray, w, 
#                                       [(0, np.radians(30), np.radians(250)),
#                                       (0, np.radians(90), np.radians(50)),
#                                       (0, np.radians(50), np.radians(145)),
#                                       (0, np.radians(50), np.radians(320))],
#                                     wavenumber)

## set plotting mode of the radiation pattern
plotMode = PlotMode.plot3d
# plotMode = PlotMode.plot2d

## range and resolution of the plot
plotThetaRange = np.pi
plotPhiRange = 2*np.pi
poleResolution = 2**6
azimuthResolution = 2**6

array_processing(antennaArray, wavenumber, w,
                poleResolution, azimuthResolution,
                plotThetaRange, plotPhiRange,
                plotMode)