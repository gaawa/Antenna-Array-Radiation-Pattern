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
from antenna_element.omnidirectional import OmnidirectionalAntennaElement
from antenna_element.halfdirectional import HalfdirectionalAntennaElement
from beamformer import no_beamformer, projection_beamformer, partial_projection_beamformer

## set physical constants
freq = 3e9
wavelength = spc.speed_of_light/freq
wavenumber = 2*np.pi/wavelength

## select antenna type
AntennaElementClass = HalfdirectionalAntennaElement
# AntennaElementClass = OmnidirectionalAntennaElement

## spawn antenna elements in linear array
antennaArray = LinearArray(AntennaElementClass, wavelength, numElements=4, elementDistanceFactor=0.5)
antennaArray = BowCircularArray(AntennaElementClass, wavelength)

## select beamformer type
beamformer = None 
beamformer = BeamFormer.Projection

## set desired beamforming direction
beamformerTheta = np.pi/5
beamformerPhi = 0

## calculate weighting factors
w = no_beamformer(antennaArray)
w = projection_beamformer(antennaArray, beamformerTheta, beamformerPhi, wavenumber)
w = partial_projection_beamformer(antennaArray, beamformerTheta, beamformerPhi, 
                                  wavenumber, active_angle=np.pi/2)

## set plotting mode of the radiation pattern
plotMode = PlotMode.plot3d

## range and resolution of the plot
plotThetaRange = np.pi
plotPhiRange = 2*np.pi
poleResolution = 2**6
azimuthResolution = 2**6

array_processing(antennaArray, wavenumber, w,
                poleResolution, azimuthResolution,
                plotThetaRange, plotPhiRange,
                plotMode)