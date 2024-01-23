import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from functools import partial
import numpy as np
import scipy.constants as spc
from antenna_pattern.circledirectional import CircledirectionalAntennaPattern
from antenna_pattern.omnidirectional import OmnidirectionalAntennaPattern
from antenna_pattern.csv_file_pattern import CsvFilePattern
from beamformer import no_beamformer, projection_beamformer, partial_projection_beamformer
from beamformer import synthesis_beamformer, synthesis_canceller_beamformer
from antenna_array import CircularArray, LinearArray
import array_processing
from plotter import plot_antenna_topology

# global variables

freq = 3e9
wavelength = spc.speed_of_light/freq
wavenumber = 2*np.pi/wavelength
nSteps = 100
nAngResolution = 200
projection = 'polar'

def animate(i, ax, bfThetas, bfPhis, antennaArrayList, legendList):
    ax.clear()
    # get current values
    bfTheta = bfThetas[i]
    bfPhi = bfPhis[i]
    
    for iArr, antennaArray in enumerate(antennaArrayList):
        # calculate beamforming weighting factors
        w = projection_beamformer(antennaArray, bfTheta, bfPhi, wavenumber)
        # w = w * np.array([[1],[0]])

        # w = synthesis_beamformer(antennaArray, [(10, bfTheta, bfPhi)], wavenumber)
        # w = len(antennaArray.arrayElements)*w/np.linalg.norm(w)
        # azimuth angle steps for simulation
        azimuthAngles = np.arange(0, 2*np.pi, 2*np.pi/nAngResolution)
        
        # calculate the 2D pattern along the azimuth angles
        pattern = array_processing.calculate_radiation_pattern(antennaArray.arrayElements, w, 
                                                               np.array([np.pi/2]), azimuthAngles,
                                                               wavenumber)
        
        patternSqueezed = np.squeeze(pattern)
        arrayPatternLog = 10*np.log10(np.abs((patternSqueezed))**2)
        if projection == 'linear':
            angAx = np.rad2deg(azimuthAngles)
        else:
            angAx = azimuthAngles
        line, = ax.plot(angAx, np.squeeze(arrayPatternLog))
        line.set_label(legendList[iArr])

    ax.set_ylim((-40, 40))
    ax.grid(True)
    ax.legend()
    ax.set_title('2 Antennas at different orientation angle')
    # ax.set_title('2 omnidirectional antennas')
    ax.set_title('2 Antennas at 135°')
    if projection == 'linear':
        xLinePos = np.rad2deg(bfPhi)
    else:
        xLinePos = bfPhi
    ax.axvline(x=xLinePos, color='r', label='steering angle')
 
if __name__ == "__main__":
    # array setup
    antennaPattern = CsvFilePattern('antenna_pattern/Gain_lin_simulation_johann.csv', debug=False, fastMode=True)
    # antennaPattern = OmnidirectionalAntennaPattern()

    antennaArrayList = []
    legendList = []
    
    # 1 antenna test 
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=1, circularAzimuthOffset=2*3/4*np.pi)
    # antennaArrayList.append(antennaArray)
    # legendList.append('0°')
    
    # 4 antennas linear array
    # antennaArray = LinearArray(antennaPattern, wavelength, numElements=4, elementDistanceFactor=0.5)
    # antennaArrayList.append(antennaArray)
    # legendList.append('linear array')

    # 4 antennas 1/2 pi arc
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=4, circularAng=2*np.pi/4, circularAzimuthOffset=np.pi/4)
    # antennaArrayList.append(antennaArray)
    # legendList.append('1/2 pi arc')

    # 4 antennas 3/4 pi arc
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=4, circularAng=3*np.pi/4, circularAzimuthOffset=np.pi/8)
    # antennaArrayList.append(antennaArray)
    # legendList.append('3/4 pi arc')

    # 2 antennas 0°
    # antennaArray = LinearArray(antennaPattern, wavelength, numElements=2)
    # antennaArrayList.append(antennaArray)
    # legendList.append('0°')

    # 2 antennas 30°
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=2, circularAng=2*np.pi/12, circularAzimuthOffset=np.pi/2-np.pi/12)
    # antennaArrayList.append(antennaArray)
    # legendList.append('30°')

    # 2 anennas 45°
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=2, circularAng=2*np.pi/8, circularAzimuthOffset=np.pi/2-np.pi/8)
    # antennaArrayList.append(antennaArray)
    # legendList.append('45°')

    # 2 antennas 60°
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=2, circularAng=2*np.pi/6, circularAzimuthOffset=np.pi/2-np.pi/6)
    # antennaArrayList.append(antennaArray)
    # legendList.append('60°')

    # 2 antennas 90°
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=2, circularAng=2*np.pi/4, circularAzimuthOffset=np.pi/2-np.pi/4)
    # antennaArrayList.append(antennaArray)
    # legendList.append('90°')

    # 2 antennas 130°
    antennaArray = CircularArray(antennaPattern, wavelength, numElements=2, circularAng=2*np.pi*3/8, circularAzimuthOffset=np.pi/2-np.pi*3/8)
    antennaArrayList.append(antennaArray)
    legendList.append('90°')

    # Linear test
    # antennaArray = LinearArray(antennaPattern, wavelength, numElements=4, elementDistanceFactor=0.5)
    # antennaArray = CircularArray(antennaPattern, wavelength, numElements=4, circularAng=0.001, circularAzimuthOffset=np.pi/2-0.0005)
    
    bfThetas = np.zeros(nSteps)
    bfThetas = np.full(nSteps, np.pi/2)
    bfPhis = np.arange(0, 2*np.pi, 2*np.pi/nSteps)
    
    # plot_antenna_topology.plot_dots(antennaArray.arrayElements)
    plot_antenna_topology.plot_quiver(antennaArray.arrayElements,
                                     ax_size=antennaArray.ax_size, 
                                     arrow_size=antennaArray.arrow_size,
                                     scale_arrow = 1)
    
    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(projection=projection)

    ani = FuncAnimation(fig,
                        partial(animate, 
                                ax=ax,
                                bfThetas=bfThetas, bfPhis=bfPhis,
                                antennaArrayList=antennaArrayList,
                                legendList=legendList),
                        frames=nSteps,
                        repeat=True)
    
    # plt.show()
    writer = PillowWriter(fps=15,
                          metadata=dict(artist='Me'),
                          bitrate=1800)
    ani.save('dual_antenna_johan_135_polar.gif', writer=writer)