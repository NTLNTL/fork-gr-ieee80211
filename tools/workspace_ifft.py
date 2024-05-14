from matplotlib import pyplot as plt

import numpy as np


def procNonDataSC(inQam):
    # input QAM has DC 0
    if (len(inQam) in [53, 117, 245]):
        return ([0] * 6 + inQam + [0] * 5)
    elif (len(inQam) == 57):
        return ([0] * 4 + inQam + [0] * 3)
    else:
        print("cloud phy80211 header, procNonDataSC: input length error %d" % (len(inQam)))
        return []
    
def procFftMod(inQam):
    # ifft shift
    if(len(inQam) in [64, 128, 256]):
        return list(np.fft.ifft((inQam[int(len(inQam)/2):] + inQam[:int(len(inQam)/2)])))
    else:
        print("cloud phy80211 header, procFftMod: input length error %d" % (len(inQam)))
        return []
    
def myPlot(inSigComp):
    plt.plot(np.real(inSigComp))
    plt.plot(np.imag(inSigComp))
    plt.show()




def myConstellationPlot(inSig):
        x = []
        y = [] 
        for i in range(len(inSig)):
            x.append(np.real(inSig))
            y.append(np.imag(inSig))
        plt.xlim([-1.1, 1.1])
        plt.ylim([-1.1, 1.1])
        plt.scatter(x,y)
        plt.show()

if __name__ == "__main__":
    # subcarrier = [0,1,1,1] + [0]*60
    # y =  np.fft.ifft(subcarrier)
    # plt.plot(np.real(y))
    # plt.show()



    C___LTF_L = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1]
    C___LTF_R = [1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1]
    C_LTF_L_26 = C___LTF_L + [0] + C___LTF_R

    print("before ifft",procNonDataSC(C_LTF_L_26))
    myConstellationPlot(procNonDataSC(C_LTF_L_26))

    ifftout = procFftMod(procNonDataSC(C_LTF_L_26))
    fftout = np.fft.fft(ifftout)
    myConstellationPlot(fftout)

    print("fft out",fftout)
    LFT = ifftout[int(64/2):] +ifftout+ifftout
    # myPlot(LFT) 

