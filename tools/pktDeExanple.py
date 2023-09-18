import os
import numpy as np
import struct
import phy80211header as p8h
from matplotlib import pyplot as plt



class dephy80211():
    def __init__(self,addr, ifaddNoise = True, ifDebug = True):
        #rx trigger 
        self.timeInSig = p8h.procLoadComplexBin(addr)
        self.threshold = 0.8
        self.autoOut = []
        noiseAmp = 0.000000000000005
        self.timeInSigAddNoise = []
        self.STFstartIndex = 0

        #adding noise 
        if ifaddNoise:
            for i in range(len(self.timeInSig)):
                real = np.real(self.timeInSig[i])
                imag = np.imag(self.timeInSig[i])
                self.timeInSigAddNoise.append(complex(real+np.random.normal() * noiseAmp, imag+np.random.normal()*noiseAmp))

        #channel info
        self.cfoRad = 0 

    def rxStepsList(self):
        #
        self.__procRxLegacyStfTrigger()

    def __procRxLegacyStfTrigger(self):
        #do auto cor 
        #return index that autoOut > threshold 
        tmpSum = 0 
        tmpSumSig1 = 0 
        tmpSumSig2 = 0
        tmpNumSearchWin = 50
        tmpWinStartIndex = 0 

        if(len(self.timeInSigAddNoise)>32):
            for i in range(16):
                tmpSum += self.timeInSigAddNoise[i]*np.conj(self.timeInSigAddNoise[i+16])
                tmpSumSig1 += np.abs(self.timeInSigAddNoise[i]) * np.abs(self.timeInSigAddNoise[i])
                tmpSumSig2 += np.abs(self.timeInSigAddNoise[i+16])* np.abs(self.timeInSigAddNoise[i+16])
        else: print("ERROR: in __procRxLegacyStfTrigger, input sig length too short")



        for i in range(len(self.timeInSigAddNoise)-32):
            self.autoOut.append(np.abs(tmpSum)/np.sqrt(tmpSumSig1)/np.sqrt(tmpSumSig2))
            tmpSum -= self.timeInSigAddNoise[i]*np.conj(self.timeInSigAddNoise[i+16])
            tmpSumSig1 -= np.abs(self.timeInSigAddNoise[i]) * np.abs(self.timeInSigAddNoise[i])
            tmpSumSig2 -= np.abs(self.timeInSigAddNoise[i+16]) * np.abs(self.timeInSigAddNoise[i+16])
            tmpSum += self.timeInSigAddNoise[i+16]*np.conj(self.timeInSigAddNoise[i+16+16])
            tmpSumSig1 += np.abs(self.timeInSigAddNoise[i+16]) * np.abs(self.timeInSigAddNoise[i+16 ])
            tmpSumSig2 += np.abs(self.timeInSigAddNoise[i+16+16]) * np.abs(self.timeInSigAddNoise[i+16+16 ])

        tmpshoulderCount = 0 
        tmpLastIndex  = 0
        for i in range(len(self.autoOut)): 
            if  self.autoOut[i]>self.threshold:
                tmpLastIndex = i
                tmpshoulderCount+= 1

        self.STFstartIndex = int(np.ceil(tmpLastIndex + tmpLastIndex-tmpshoulderCount)/2) - 63
        print("self.STFstartIndex",self.STFstartIndex)
        myPlot(self.timeInSigAddNoise)



def myPlot(inSigComp):
    plt.plot(np.real(inSigComp))
    plt.plot(np.imag(inSigComp))
    plt.show()

if __name__ == "__main__":



    pyToolPath = os.path.dirname(__file__)
    addr = os.path.join(pyToolPath, "../tmp/sig80211GenMultipleSiso_1x1_0.bin")
    p = dephy80211(addr,True,False)
    p.rxStepsList()



