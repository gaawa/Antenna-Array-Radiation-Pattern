from antenna_pattern.antenna_pattern_abs import AntennaPattern
import numpy as np
import pandas
import utils

import matplotlib.pyplot as plt

class CsvFilePattern(AntennaPattern):
    """
    Reads antenna radiation pattern from a file. Values are converted to linear scale for output

    TODO: angle offset parameter to align the pattern
    TODO: interpolate inbetween angles
    """
    def __init__(self, filePath, debug=False, fastMode=True):
        super().__init__()

        # csv data frame
        self.radiationPatternValues = pandas.read_csv(filePath)
        self.thetaCsvAxis = np.asarray(self.radiationPatternValues["Theta[deg]"])*np.pi/180
        self.thetaCsvAxisUnique = np.unique(np.asarray(self.radiationPatternValues["Theta[deg]"]))*np.pi/180
        self.phiCsvAxis = np.asarray(self.radiationPatternValues["Phi[deg]"])*np.pi/180
        self.phiCsvAxisUnique = np.unique(np.asarray(self.radiationPatternValues["Phi[deg]"]))*np.pi/180
        try:
            self.gain = np.asarray(self.radiationPatternValues["GainTotal"])
        except KeyError:
            self.gain = 10**(np.asarray(self.radiationPatternValues["dB(GainTotal)"])/10)
        # Reshape gain array to matrix form
        # axis 0: phi axis
        # axis 1: theta axis
        self.gainMatrix = self.gain.reshape((self.phiCsvAxisUnique.size, self.thetaCsvAxisUnique.size), order='F')
        
        self.fastMode = fastMode
        self.debug = debug

   
    def get_element_factor(self, signalVecLocal):
        if self.fastMode:
            theta = signalVecLocal[0,:]
            phi = signalVecLocal[1,:]
            
            # get closest array value indexes for each angles
            idxTheta = utils.get_closest_avi(self.thetaCsvAxisUnique, theta)
            idxPhi = utils.get_closest_avi(self.phiCsvAxisUnique, phi)
            
            pattern = self.gainMatrix[idxPhi, idxTheta]
        else:
            pattern = np.zeros(len(signalVecLocal[0,:]))
            for pos in range(len(signalVecLocal[0,:])):

                error = (np.array([self.thetaCsvAxis, self.phiCsvAxis])
                         - signalVecLocal[:, pos][:,None].astype(float) )**2
                #error = (np.exp(1j*theta)*np.conj(1j*signalVecLocal[0, pos]) + np.exp(1j*phi)*np.conj(1j*signalVecLocal[1, pos])).real

                patternIndex = np.argmin(np.sum(error, axis=0))
                pattern[pos] = self.gain[patternIndex]
                #print(patternIndex)

        if self.debug:
            fig, ax = plt.subplots(1,2)
            ax[0].plot(signalVecLocal[0,:])
            ax[0].plot(signalVecLocal[1,:])
            ax[1].plot(self.thetaCsvAxis)
            ax[1].plot(self.phiCsvAxis)
            ax[0].legend(["theta","phi"])
            ax[1].legend(["theta","phi"])
            ax[0].set_title(["Program Winkel"])
            ax[1].set_title(["CSV Winkel"])

            plt.show(block=True)

            N = np.int_(np.sqrt(signalVecLocal.shape[1])) # N x N antenna array pattern 

            pattern = pattern.reshape((N,N))

            plt.imshow(pattern)
            plt.show(block=True)
            
            plt.imshow(self.gainMatrix)
            plt.show(block=True)

            pattern = pattern.reshape(N*N)

        return pattern