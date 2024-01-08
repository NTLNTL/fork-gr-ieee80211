from matplotlib import pyplot as plt

import numpy as np
if __name__ == "__main__":
    subcarrier = [0,1,1,1] + [0]*60
    y =  np.fft.ifft(subcarrier)
    plt.plot(np.real(y))
    plt.show()
    