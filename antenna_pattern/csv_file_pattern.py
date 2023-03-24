from antenna_pattern.antenna_pattern_abs import AntennaPattern
import numpy as np
import pandas

import matplotlib.pyplot as plt

class CsvFilePattern(AntennaPattern):
    """
    Reads antenna radiation pattern from a file. Values are converted to linear scale for output
    """
    def __init__(self, filePath, debug=False):
        super().__init__()
        self.radiationPatternValues = pandas.read_csv(filePath)
        self.debug = debug

   
    def get_element_factor(self, signalVecLocal):
        pattern = np.zeros(len(signalVecLocal[0,:]))

        # Theta index in radians
        theta = np.asarray(self.radiationPatternValues["Theta[deg]"]*np.pi/180)
        # Phi index in radians      
        phi = np.asarray(self.radiationPatternValues["Phi[deg]"]*np.pi/180) 
        # Linear antenna gain
        gain = self.radiationPatternValues["GainTotal"]


        for pos in range(len(signalVecLocal[0,:])):

            error = (np.array([theta, phi])-signalVecLocal[:, pos][:,None].astype(float))**2
            #error = (np.exp(1j*theta)*np.conj(1j*signalVecLocal[0, pos]) + np.exp(1j*phi)*np.conj(1j*signalVecLocal[1, pos])).real

            patternIndex = np.argmin(np.sum(error, axis=0))
            pattern[pos] = gain[patternIndex]
            #print(patternIndex)

        if self.debug:
            fig, ax = plt.subplots(1,2)
            ax[0].plot(signalVecLocal[0,:])
            ax[0].plot(signalVecLocal[1,:])
            ax[1].plot(theta)
            ax[1].plot(phi)
            ax[0].legend(["theta","phi"])
            ax[1].legend(["theta","phi"])
            ax[0].set_title(["Program Winkel"])
            ax[1].set_title(["CSV Winkel"])

            plt.show(block=True)

            N = np.int_(np.sqrt(signalVecLocal.shape[1])) # N x N antenna array pattern 

            pattern = pattern.reshape((N,N))

            plt.imshow(pattern)
            plt.show(block=True)

            pattern = pattern.reshape(N*N)

        return pattern