import os
import numpy as np
import math
import phy80211header as p8h
from matplotlib import pyplot as plt
import phy80211 



class dephy80211():
    def __init__(self,addr, ifaddNoise = True, ifDebug = True):
        #import
        dephy80211 = phy80211.phy80211(ifDebug=False)

        #rx trigger 
        self.timeInSig = p8h.procLoadComplexBin(addr)
        self.threshold = 0.8
        self.searchWin = 50
        self.autoOut = []
        noiseAmp = 0.0001
        self.STFstartIndex = 0


        #adding noise prevent denominator is zero
        if ifaddNoise:
            for i in range(len(self.timeInSig)):
                real = np.real(self.timeInSig[i])
                imag = np.imag(self.timeInSig[i])
                self.timeInSig[i] = complex(real+np.random.normal() * noiseAmp, imag+np.random.normal()*noiseAmp)
        #channel info
        self.cfoRad = 0 


        # debug
        self.ifdb = ifDebug

    def rxStepsList(self):
        self.__procRxLegacyStfTrigger()
        self.__procRxLegacyCfoEst()

    def __procRxLegacyStfTrigger(self):
        #do auto cor 
        #return index that autoOut > threshold 
        tmpSum = 0 
        tmpSumSig1 = 0 
        tmpSumSig2 = 0
        if(len(self.timeInSig)>32):
            for i in range(16):
                tmpSum += self.timeInSig[i]*np.conj(self.timeInSig[i+16])
                tmpSumSig1 += np.abs(self.timeInSig[i]) * np.abs(self.timeInSig[i])
                tmpSumSig2 += np.abs(self.timeInSig[i+16])* np.abs(self.timeInSig[i+16])
        else: print("ERROR: in __procRxLegacyStfTrigger, input sig length too short")
    





        for i in range(len(self.timeInSig)-32):
            if np.real(tmpSumSig1)!= 0 and np.real(tmpSumSig2)!= 0:
                self.autoOut.append(np.abs(tmpSum)/np.sqrt(tmpSumSig1)/np.sqrt(tmpSumSig2))
            else:
                self.autoOut.append(0)
            tmpSum -= self.timeInSig[i]*np.conj(self.timeInSig[i+16])
            tmpSumSig1 -= np.abs(self.timeInSig[i]) * np.abs(self.timeInSig[i])
            tmpSumSig2 -= np.abs(self.timeInSig[i+16]) * np.abs(self.timeInSig[i+16])
            tmpSum += self.timeInSig[i+16]*np.conj(self.timeInSig[i+16+16])
            tmpSumSig1 += np.abs(self.timeInSig[i+16]) * np.abs(self.timeInSig[i+16 ])
            tmpSumSig2 += np.abs(self.timeInSig[i+16+16]) * np.abs(self.timeInSig[i+16+16 ])
           
        tmpshoulderCount = 0 
        tmpLastIndex  = 0
        for i in range(len(self.autoOut)): 
            if  self.autoOut[i]>self.threshold:
                tmpLastIndex = i
                tmpshoulderCount+= 1

        self.STFstartIndex = int(np.ceil(tmpLastIndex + tmpLastIndex-tmpshoulderCount)/2) - 63
        if self.ifdb:print("self.STFstartIndex",self.STFstartIndex)

        
        
    def __procRxLegacyCfoEst(self):
        inSig = self.timeInSig[self.STFstartIndex:self.STFstartIndex+320]

        #estimate tmpCoarseCfo from STF
        tmpCoarseCfo = np.mean([inSig[i] * np.conj(inSig[i + 16]) for i in range(0, 160-16)])
        tmpCoarseCfoHz = np.arctan2(np.imag(tmpCoarseCfo), np.real(tmpCoarseCfo)) / 16 * 20000000 / 2 / np.pi
        tmpCoarseCfoRad = np.arctan2(np.imag(tmpCoarseCfo), np.real(tmpCoarseCfo)) / 16

        #compensate cfo in LTF by using tmpCoarseCfo
        tmpLTF =  []
        for i in range (160):
            tmpLTF.append(inSig[i+160] * complex(np.cos(i *tmpCoarseCfoRad ), np.sin(i * tmpCoarseCfoRad)))
        # myConstellationPlot(np.fft.fft(tmpLTF))
        tmpFineCfo = np.mean([tmpLTF[i] * np.conj(tmpLTF[i + 64]) for i in range(0, 64)])
        tmpFineCfoRad = np.arctan2(np.imag(tmpFineCfo), np.real(tmpFineCfo))/64

        self.cfoRad = tmpCoarseCfoRad + tmpFineCfoRad
        if self.ifdb:print("self.__procRxLegacyCfoEst:",self.cfoRad* 20000000 / 2 / np.pi)

            
                                                                                 











def myConstellationPlot(inSig):
        x = []
        y = [] 
        for i in range(len(inSig)):
            x.append(np.real(inSig))
            y.append(np.imag(inSig))
        # plt.xlim([-1.1, 1.1])
        # plt.ylim([-1.1, 1.1])
        plt.scatter(x,y)
        plt.show()



def myPlot(inSigComp):
    plt.plot(np.real(inSigComp))
    plt.plot(np.imag(inSigComp))


if __name__ == "__main__":



    pyToolPath = os.path.dirname(__file__)
    addr = os.path.join(pyToolPath, "../tmp/sig80211GenMultipleSiso_1x1_0.bin")
    # addr =   "C:\\Users\\naton\\Dropbox\\PhD\\80211\\code_80211a\sig80211LegacyGenOne.bin_1x1_0.bin"
    p = dephy80211(addr,ifaddNoise = False,ifDebug = True )
    p.rxStepsList()



