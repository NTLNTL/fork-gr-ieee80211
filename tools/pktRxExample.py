import os
import numpy as np
import math
import phy80211header as p8h
import mac80211header as m8h
from matplotlib import pyplot as plt
import cmath
import struct
import mac80211 as mac
import txbits
import inspect
class dephy80211siso():
    def __init__(self,addr, ifaddNoise = True, ifDebug = True):


        #rx trigger 
        self.timeInSig = p8h.procLoadComplexBin(addr)
        
        self.threshold = 0.8
        self.searchWin = 50
        self.autoOut = []
        noiseAmp = 0.0001
        self.legacyStfIndex = 0

        #adding noise prevent denominator is zero
        if ifaddNoise:
            for i in range(len(self.timeInSig)):
                real = np.real(self.timeInSig[i])
                imag = np.imag(self.timeInSig[i])
                self.timeInSig[i] = complex(real+np.random.normal() * noiseAmp, imag+np.random.normal()*noiseAmp)
       
        #channel info
        self.cfoRad = 0 
        self.chanInfoLtf = []

        #rx legacy sig field info 
        self.LSigBits = []
        self.mdpuLen = 0

        #rx vhtSigB field 
        self.vhtSigBBits = None
        self.calVhtSigBCrc = 0 

        #rx ht sig field 
        self.aggre = 0


        #create rx demod class from received pkta
        self.demod = None
        self.nSym = 0

        #rx data sym
        self.pilotPIdx = 0  # P_0 is for signal symbol, P_1 starts from DATA symbol (legacy);
        self.nSymProcIdx = 0 #include both L-STF, L-LTF, L-SIG (5 nSym).   
        self.rxDataSym = None
        self.rxDataQam = []
        self.rxInterleaveBits = []
        self.rxCodedBits = []
        self.rxScrambleBits = []
        self.rxDatabits = []
        self.nSS = 1 
        self.mcs = 0 
        self.phyFormat = None

        #rx phy bits out 
        self.phyBits = []
        # debug
        self.ifdb = ifDebug

    def rxSisoStepsList(self): 
        # self.__testNoCfo()
        self.__procRxLegacyStfTrigger()
        self.__procRxLegacyCfoEst()
        self.__procRxLegacyChanEst()
        self.__procRxLegacySigSym()
        self.__procRxPktFormat()
        if self.nSS == 1: # siso 
            if self.phyFormat == p8h.F.L:
                self.__procRxSisoData()
            elif self.phyFormat == p8h.F.VHT:
                self.__procRxChanUpdate()
                self.__procRxParseVHTSigB()
                self.__procRxSisoData()
            elif self.phyFormat == p8h.F.HT:
                self.__procRxChanUpdate()
                self.__procRxSisoData()
        else:
            print("in rxSisoStepsList self.nSS  ",self.nSS) 
            print("__procRxPassSymToMimoRx len",len(self.__procRxPassSymToMimoRx()) ) #should all go through cfo, symbols start from VHT-STF,VHT-LTF has new channel estimation
            return self.__procRxPassSymToMimoRx() #return time domain 64 samples for each symbols
 
        return self.rxDatabits
    

    def __procRxParseVHTSigB(self):
        
        tmpSigBMod = p8h.modulation(phyFormat=self.phyFormat, mcs=0, bw=p8h.BW.BW20, nSTS=1, shortGi=False) #mcs set to 0 for VHTSigB interleave
        tmpRxVhtBSym = self.timeInSig[self.legacyStfIndex+(self.nSymProcIdx*80):self.legacyStfIndex+(self.nSymProcIdx*80)+80]
        tmp = []
        for nSampleIter in range(80):
            tmp.append(tmpRxVhtBSym[nSampleIter]*complex(np.cos(((self.nSymProcIdx*80)+nSampleIter)*self.cfoRad), np.sin(((self.nSymProcIdx*80)+nSampleIter)*self.cfoRad)))
        # remove CP 
        tmpSymTime = tmp[16:80]
        tmpSigFreq = p8h.procFftDemod(tmpSymTime)
        tmpSigFreq = self.procCompChan(tmpSigFreq) #64
        tmpSigFreq = p8h.procRmDcNonDataSc(tmpSigFreq, self.phyFormat)
        tmpSigFreq = self.procPilotTrack(tmpSigFreq,p8h.C_PILOT_VHT[p8h.BW.BW20.value])
        tmpLlr = p8h.procRemovePilots(tmpSigFreq) 
        tmpSigBMod.nSym = 1 # for interleave
        tmpCoded = self.procDeinterleave(list(np.real(tmpLlr)),tmpSigBMod)
        tmpVhtSigBBits = p8h.procViterbiDecoder(tmpCoded,26, p8h.CR.CR12)
        self.calVhtSigBCrc = p8h.genBitBitCrc8(tmpVhtSigBBits[0:20]) #crc bits for vhtSigB calculated but not verify in service 
        self.vhtSigBBits = tmpVhtSigBBits 
        self.mdpuLen = self.bi2Deci(self.vhtSigBBits[0:17]) * 4
        if self.ifdb: print("DEBUG:VHT sigB bits", tmpVhtSigBBits,self.calVhtSigBCrc)
        
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
        # myPlot(self.autoOut)
        self.legacyStfIndex = int(np.ceil(tmpLastIndex + tmpLastIndex-tmpshoulderCount)/2) - 63
        self.legacyStfIndex = 1200
        if self.ifdb:print("DEBUG: self.legacyStfIndex",self.legacyStfIndex)
  
    def __procRxLegacyCfoEst(self):
        inSig = self.timeInSig[self.legacyStfIndex:self.legacyStfIndex+320]

        #estimate tmpCoarseCfo from STF
        tmpCoarseCfo = np.mean([inSig[i] * np.conj(inSig[i + 16]) for i in range(0, 160-16)])
        tmpCoarseCfoHz = np.arctan2(np.imag(tmpCoarseCfo), np.real(tmpCoarseCfo)) / 16 * 20000000 / 2 / np.pi
        tmpCoarseCfoRad = np.arctan2(np.imag(tmpCoarseCfo), np.real(tmpCoarseCfo)) / 16

        #compensate cfo in LTF by using tmpCoarseCfo
        tmpLTF =  []
        for i in range (160):
            tmpLTF.append(inSig[i+160] * complex(np.cos(i *tmpCoarseCfoRad), np.sin(i * tmpCoarseCfoRad ))) # check 60 and 61 sample offset

        
        tmpFineCfo = np.mean([tmpLTF[i] * np.conj(tmpLTF[i + 64]) for i in range(0, 160-64)])
        tmpFineCfoRad = np.arctan2(np.imag(tmpFineCfo), np.real(tmpFineCfo))/64

        #compensate cfo in LTF by using tmpFineCfoRad
        for i in range (160): 
            tmpLTF[i] = (tmpLTF[i] * complex(np.cos(i *tmpFineCfoRad ), np.sin(i * tmpFineCfoRad  )))
        self.cfoRad = tmpCoarseCfoRad + tmpFineCfoRad
        self.timeInSig[self.legacyStfIndex+160:self.legacyStfIndex+160+160] = tmpLTF




        self.nSymProcIdx = 4
        if self.ifdb:print("DEBUG: self.cfoRad in frequency:",self.cfoRad,self.cfoRad* 20000000 / 2 / np.pi)

    def __procRxLegacyChanEst(self):# 64 samples channel info 
                                    #channel before ifft compensate, 
                                    # after channel compensate, process symbol
        tmpChan = []
        tmpLtfOri = p8h.procNonDataSC(p8h.C_LTF_L[p8h.BW.BW20.value])
        #flip by half 
        tmpRxLtf1 = p8h.procFftDemod(self.timeInSig[self.legacyStfIndex+160+32:self.legacyStfIndex+160+32+64])
        tmpRxLtf2 = p8h.procFftDemod(self.timeInSig[self.legacyStfIndex+160+32+64:self.legacyStfIndex+160+32+64+64])


        for i in range(0,len(tmpRxLtf1)):
            if(tmpLtfOri[i] == 0):
                tmpChan.append(0)
            else:
                tmpChan.append((tmpRxLtf1[i]+tmpRxLtf2[i])/2/tmpLtfOri[i])
                # tmpChan.append((tmpRxLtf1[i])/tmpLtfOri[i])
        self.chanInfoLtf = tmpChan
        # compensate channel
        # myPlot(np.fft.ifftself.procCompChan(tmpRxLtf1))
        
        # self.procCompChan(tmpRxLtf1)
        # tmpOut = []
        # for i in range(0,len(tmpRxLtf1)):
        #     if(tmpChan[i] == 0):
        #         tmpOut.append(0)
        #     else:
        #         tmpOut.append(tmpRxLtf1[i]/tmpChan[i])
        # myConstellationPlot(tmpOut)
    def __procRxChanUpdate(self):
        self.nSymProcIdx = 8
        tmpLTF = self.timeInSig[self.legacyStfIndex+(self.nSymProcIdx*80):self.legacyStfIndex+(self.nSymProcIdx*80)+80]
        tmp = []
        tmpLtfOri = p8h.procNonDataSC(p8h.C_LTF_HT[p8h.BW.BW20.value])
        for nSampleIter in range(80):
            tmp.append(tmpLTF[nSampleIter]*complex(np.cos((self.nSymProcIdx * 80 + nSampleIter)*self.cfoRad), np.sin((self.nSymProcIdx * 80 + nSampleIter)*self.cfoRad)))
        # remove CP 
        tmpSymTime = tmp[16:80]
        tmpSymFreq = p8h.procFftDemod(tmpSymTime)
        tmpChan = [0 if tmpLtfOri[i] == 0 else tmpSymFreq[i]/tmpLtfOri[i] for i in range(len(tmpLtfOri)) ]
        self.chanInfoLtf = tmpChan
        self.nSymProcIdx += 1
        print("channel updated")

    def __procRxSisoData(self):



        if self.phyFormat == p8h.F.L:
            self.nSymProcIdx = 5
        elif  self.phyFormat == p8h.F.HT:
            self.nSymProcIdx = 9 
        elif self.phyFormat == p8h.F.VHT:
            self.nSymProcIdx = 10

        self.rxDataSym = self.timeInSig[self.legacyStfIndex+(self.nSymProcIdx*80):self.legacyStfIndex+(self.nSymProcIdx*80)+self.demod.nSym*80]

        # myPlot(self.rxDataSym)
        tmpPilot = p8h.C_PILOT_L
        for nSymIter in range(self.demod.nSym):
            tmpSym = self.rxDataSym[nSymIter*80:nSymIter*80+80]
            tmp = []
            # compensate cfo 
            for nSampleIter in range(80):
                tmp.append(tmpSym[nSampleIter]*complex(np.cos(((nSymIter*80)+nSampleIter)*self.cfoRad), np.sin(((nSymIter*80)+nSampleIter)*self.cfoRad)))
            # remove CP 
                tmpSymTime = tmp[16:80]

            #convert to freq domain and recover channel
            tmpSigFreq = p8h.procFftDemod(tmpSymTime)
            tmpSigFreq = self.procCompChan(tmpSigFreq) #64
            tmpSigFreq = p8h.procRmDcNonDataSc(tmpSigFreq, self.demod.phyFormat)
            tmpSigFreq = self.procPilotTrack(tmpSigFreq,tmpPilot)
            myConstellationPlot(tmpSigFreq)
            print("len procPilotTrack", len(tmpSigFreq))
            self.rxDataQam += p8h.procRemovePilots(tmpSigFreq) 
            # myConstellationPlot(p8h.procRemovePilots(tmpSigFreq) )
            print("len  self.rxDataQam", len(self.rxDataQam))
            self.nSymProcIdx +=1 
            if(not(self.phyFormat == p8h.F.L)):
                tmpPilot = tmpPilot[1:] + [tmpPilot[0]] #?
            # myConstellationPlot( self.rxDataQam)



        if self.ifdb:print("DEBUG: self.rxDataQam len:",len(self.rxDataQam))
        # if self.ifdb:print("DEBUG: self.rxDataQam", self.rxDataQam)
        self.rxInterleaveBits = self.procSymQamToLlr()
        # print("compare interleave",compare([1 if np.real(each)>0 else 0  for each in self.rxInterleaveBits ],txbits.htMcs0interBits))
        # if self.ifdb:print("DEBUG: self.rxInterleaveBits len:",len(self.rxInterleaveBits))
        # if self.ifdb:print("DEBUG: self.rxInterleaveBits: \n", self.rxInterleaveBits)

        self.rxCodedBits = self.procDeinterleave(self.rxInterleaveBits,self.demod)
        if self.ifdb:print("DEBUG: self.rxCodedBits len:",len(self.rxCodedBits))
        # print("compare coded",compare( [1 if np.real(each)>0 else 0  for each in self.rxCodedBits ],txbits.htMcs0CodedBits))
        # if self.ifdb:print("DEBUG: self.rxCodedBits: \n", self.rxCodedBits)
        
        self.rxScrambleBits = p8h.procViterbiDecoder(self.rxCodedBits,self.demod.nDBPS*self.demod.nSym,self.demod.cr)

        if self.ifdb:print("DEBUG: self.rxScrambleBits len:",len(self.rxScrambleBits))
        # if self.ifdb:print( "DEBUG: self.rxScrambleBits",self.rxScrambleBits)

        self.rxDatabits = self.PrcoDeScrambleBits()
        print("data bits len",len(self.rxDatabits))
        # print("data compare", compare(self.rxDatabits,txbits.vhtDataBitsMcs7))



        if self.ifdb:print("DEBUG: self.rxDatabits len:", len(self.rxDatabits),self.rxDatabits[0:16]   )
        # if self.ifdb:print("DEBUG: self.rxDatabits \n", self.rxDatabits )
    def __procRxPassSymToMimoRx(self): # return HT-STF, HT-LTF(s),Data(s) and  apply cfo 
        self.nSymProcIdx = 7
        tmpnSym = 1 + self.nSS + self.demod.nSym  #HT-STF+HT-LTF(s)+Data(s)
        tmpRxMimoSym = self.timeInSig[self.legacyStfIndex+(self.nSymProcIdx*80):self.legacyStfIndex+(self.nSymProcIdx*80)+tmpnSym*80]
        print("in__procRxPassSymToMimoRx",self.demod.nSym,tmpnSym,self.nSS)

        tmpOut = []
        for nSymIter in range(tmpnSym):
            # compensate cfo 
            tmpSym =tmpRxMimoSym[nSymIter*80:nSymIter*80+80]
            tmp = []
            for nSampleIter in range(80):
                tmp.append(tmpSym[nSampleIter]*complex(np.cos(((nSymIter*80)+nSampleIter)*self.cfoRad), np.sin(((nSymIter*80)+nSampleIter)*self.cfoRad)))
                # remove CP 
            tmpSymTime = tmp[16:80]
            tmpOut+=tmpSymTime
        return tmpOut

                



        
   
    def PrcoDeScrambleBits(self):
        tmpScrambler = 0
        tmpDeScrambledBits = []

        for i in range(7):
            tmpDeScrambledBits.append(0)
            if self.rxScrambleBits[i] == 1:
                tmpScrambler |= 1 << (6-i)
        tmpScrambler = tmpScrambler
        for i in range(7,len(self.rxScrambleBits) ):
            tmpFeedback = int(not not(tmpScrambler & 64)) ^ int(not not(tmpScrambler & 8))

            tmpDeScrambledBits.append(tmpFeedback ^ self.rxScrambleBits[i])
            tmpScrambler = ((tmpScrambler << 1) & 0x7e) | tmpFeedback  #0x7e - 126
        print("self.mdpuLen",self.mdpuLen)
        # return tmpDeScrambledBits[16:self.mdpuLen * 8 +16]#removing 16 serving bits and 6 tail bits  

        if self.phyFormat == p8h.F.VHT:
            if tmpDeScrambledBits[8:16] == self.calVhtSigBCrc:
                print("vht SigB crc check pass")
                tmpAmpduBits = tmpDeScrambledBits[16:16+8*self.mdpuLen]
                C_VHT_EOF =  tmpDeScrambledBits[16+8*self.mdpuLen:16+8*self.mdpuLen+len(p8h.C_VHT_EOF * self.demod.nPadEof)]
                print(len(C_VHT_EOF), print(C_VHT_EOF[0:len(p8h.C_VHT_EOF )]))
                return tmpAmpduBits
            else:
                print("ERROR: vht SigB crc check fail")
                tmpDeScrambledBits = []
        else: 
            tmpDeScrambledBits = tmpDeScrambledBits[16:self.mdpuLen * 8 +16]#removing 16 serving bits and 6 tail bits  
            return tmpDeScrambledBits
        # legacy and ht: data bits + tail + padding -> scramble
        # vht: vht gets data bits scrambled first then add tail bits and then coded
    def procDeinterleave(self,interleavedBits,mod):
        
        if self.phyFormat == p8h.F.L:
            tmpDeinterleaveBits = [0] * mod.nSym * mod.nCBPS
            print("in procDeinterleave",mod.phyFormat, mod.cr, mod.nSym * mod.nCBPS, mod.nCBPS)
            s = int(max(1, mod.nBPSCS/2))
            for symIter in range(mod.nSym):
                for j in range(mod.nCBPS):
                    i = int( int(s * np.floor(j/s)) +  int((j+np.floor(16*j/mod.nCBPS))%s))  
                    k = int(16 * i - (mod.nCBPS - 1) * int(np.floor(16 * i /mod.nCBPS)))
                    tmpDeinterleaveBits[symIter*mod.nCBPS+k] =  interleavedBits[symIter*mod.nCBPS+j]
            
        elif self.phyFormat== p8h.F.HT or self.phyFormat== p8h.F.VHT:
            print("in procDeinterleave ht ",mod.nSym, mod.cr,mod.nCBPSS)
            tmpDeinterleaveBits =  [0] * mod.nSym * mod.nCBPSS
            s = int(max(1, mod.nBPSCS/2))
            for symIter in range(mod.nSym):#mod.nSym
                if self.nSS == 1:
                    for r in range(mod.nCBPSS):
                        j = r
                        i = s * int(np.floor(j/s)) + (j+ int(np.floor(mod.nIntlevCol*j/mod.nCBPSS)))%s
                        k = mod.nIntlevCol * i - (mod.nCBPSS - 1 ) * int( np.floor(i / mod.nIntlevRow))
                        tmpDeinterleaveBits[symIter*mod.nCBPSS+k] = interleavedBits[symIter*mod.nCBPSS+r]
        return tmpDeinterleaveBits  
    def procSymQamToLlr(self):
        tmpBiOut = []
        tmpKmod = 0
        
        tmpDeMapping = [0] * (self.demod.nSym * self.demod.nCBPS)
        
        print("in procSymQamToLlr",self.demod.nSym * self.demod.nCBPS, self.demod.nSD)
        if self.demod.mod == p8h.M.BPSK:
            tmpKmod = 1
            for i in range (self.demod.nSym * self.demod.nSD):
                tmp = self.rxDataQam[i]
                tmpDeMapping[i] = np.real(tmp)

        elif self.demod.mod == p8h.M.QBPSK:
            tmpKmod = 1
            for i in range (self.demod.nSym * self.demod.nSD):
                tmp = self.rxDataQam[i]
                tmpDeMapping[i] = np.imag(tmp) 


        elif self.demod.mod ==  p8h.M.QPSK:
            
            tmpKmod = math.sqrt(2)
            for i in range (self.demod.nSym * self.demod.nSD):
                tmp = self.rxDataQam[i] * tmpKmod
                tmpDeMapping[i*2] = np.real(tmp)
                tmpDeMapping[i*2+1] = np.imag(tmp)
            

        elif self.demod.mod ==  p8h.M.QAM16:
            tmpKmod = math.sqrt(10)
            for i in range (self.demod.nSym * self.demod.nSD):
                tmp = self.rxDataQam[i] * tmpKmod
                tmpDeMapping[i*4] = np.real(tmp)
                tmpDeMapping[i*4+1] = 2 - np.abs(np.real(tmp))
                tmpDeMapping[i*4+2] = np.imag(tmp)
                tmpDeMapping[i*4+3] = 2 - np.abs(np.imag(tmp))


        elif self.demod.mod == p8h.M.QAM64:
            tmpKmod = math.sqrt(42)
            for i in range (self.demod.nSym * self.demod.nSD):
                tmp = self.rxDataQam[i] * tmpKmod
                tmpDeMapping[i*6] = np.real(tmp)
                tmpDeMapping[i*6+1] = 4 - abs(np.real(tmp))
                tmpDeMapping[i*6+2] = 2 - abs (4 - abs(np.real(tmp)))
                tmpDeMapping[i*6+3] = np.imag(tmp)
                tmpDeMapping[i*6+4] = 4 - abs(np.imag(tmp))
                tmpDeMapping[i*6+5] = 2 - abs (4 - abs(np.imag(tmp)))
            
        elif self.demod.mod == p8h.M.QAM256:
            pass
        elif self.demod.mod == p8h.M.QAM1024:
            pass

        # for  i in range(len(tmpDeMapping)): # for debugung, binary output
        #     if(tmpDeMapping[i] > 0):
        #         tmpBiOut.append(1)
        #     else:
        #         tmpBiOut.append(0)
        
        # print("lnterleaveOutBi", tmpBiOut)
        # print("checking", compare(tmpBiOut,tx_inter_5))
        return tmpDeMapping
    
    def __procRxPktFormat(self):
        tmpSig1 = self.timeInSig[self.legacyStfIndex+5*80:self.legacyStfIndex+5*80+80]
        tmpSig2 = self.timeInSig[self.legacyStfIndex+6*80:self.legacyStfIndex+6*80+80]
        if(sum(self.LSigBits[0:17])%2 == self.LSigBits[17] ): #parity check for Sig field
            print("pass parity check")
            if self.mcs == 0:
                print("pilot num for HTsig",self.pilotPIdx)
                self.pilotPIdx = 1 
                self.nSymProcIdx = 4  
                tmpHtSigBits = p8h.procViterbiDecoder(self.proRxSigSym(tmpSig1,p8h.M.QBPSK) + self.proRxSigSym(tmpSig2,p8h.M.QBPSK),48,  p8h.CR.CR12)
                self.pilotPIdx = 1 
                self.nSymProcIdx = 4
                tmpVhtSigBits = p8h.procViterbiDecoder(self.proRxSigSym(tmpSig1,p8h.M.BPSK) + self.proRxSigSym(tmpSig2,p8h.M.QBPSK),48,  p8h.CR.CR12)
                if (self.htVhtCrc8Check(tmpHtSigBits)):
                    if self.ifdb:print("DEBUG: HT-SIG bits:",tmpHtSigBits)
                    self.pilotPIdx = 3
                    print("ht sig fieldhtCrc8Check pass")
                    print("len HT", self.bi2Deci(tmpHtSigBits[8:24]))
                    self.mdpuLen = self.bi2Deci(tmpHtSigBits[8:24])
                    self.phyFormat = p8h.F.HT 
                    self.mcs = self.bi2Deci(tmpHtSigBits[0:7])
                    self.nSS = int(np.floor(self.mcs / 8)) + 1 #ht has 8-16 to identify 2 nSS, # vht?

                if (self.htVhtCrc8Check(tmpVhtSigBits)):
                    if self.ifdb:print("DEBUG: VHT-SIGA  bits:",tmpHtSigBits)
                    print("vht sigA field vhtCrc8Check pass")
                    self.pilotPIdx = 3
                    self.phyFormat = p8h.F.VHT
                    self.mcs = self.bi2Deci(tmpVhtSigBits[28:32])
                    self.nSS = p8h.N_STS_VHTSU[self.bi2Deci(tmpVhtSigBits[10:13])]

                if (not (self.htVhtCrc8Check(tmpHtSigBits) or self.htVhtCrc8Check(tmpVhtSigBits))):
                    print("state machin to L mcs 0 ")
                    self.pilotPIdx = 1 
                    self.phyFormat = p8h.F.L
            else: #demodulate legacy mcs > 0 
                self.phyFormat = p8h.F.L
                self.pilotPIdx = 1 
                
            print("self.nSS in pktformat",  self.nSS )
            if self.ifdb:print("DEBUG: final format:",self.phyFormat)
            if self.ifdb:print("DEBUG: final mcs:",self.mcs) 
            self.demod = p8h.modulation(self.phyFormat, self.mcs, bw=p8h.BW.BW20, nSTS=self.nSS, shortGi=False)

            if self.phyFormat ==  p8h.F.L or self.phyFormat ==  p8h.F.HT: 
                self.demod.procPktLenNonAggre(self.mdpuLen)
            else:
                self.demod.procPktLenAggre(self.mdpuLen)
                print("__procRxPktFormat VHT", p8h.F.VHT ,self.mdpuLen,self.demod.nSym )
                # self.nSymProcIdx start of data symbol 

        else:
            return "ERROR: legacy Sig parity check fail"
            
    def procPilotTrack(self,inQam,IN_PILOT): #procResiCfoCompensate 
        tmpOut = None
        tmpPilotLocation = None
        tmpStandardPilot = [int (each * p8h.C_PILOT_PS[self.pilotPIdx]) for each in IN_PILOT ]
        print("IN_PILOT",p8h.C_PILOT_PS[self.pilotPIdx],tmpStandardPilot)
        # print("pilot number in procPilotTrack", self.pilotPIdx, p8h.C_PILOT_PS[self.pilotPIdx],tmpStandardPilot)
        if len(inQam) == 52: 
            tmpPilotLocation= [5,19,32,46]
            # print("receive pilot",[inQam[each] for each in tmpPilotLocation] )
            
        elif len(inQam) == 56: 
            tmpPilotLocation = [7,21,34,48]
            # print("receive pilot",[inQam[each] for each in tmpPilotLocation] )

        elif len(inQam) == 116: #40BW, need verify 
            tmpPilotLocation = [5,21,47,68,82,110]
    
        elif len(inQam) == 244: #80BW,need verify 
            tmpPilotLocation = [19,47,83,111,132,160,196,224]

        else:
            
            return ("ERROR: procPilotTrack len of inQam error ")
            
        tmpDiff = 0
        for i in range(len(tmpPilotLocation)):
            tmpDiff += (inQam[tmpPilotLocation[i]] * np.conj( tmpStandardPilot[i]))         
        tmpAbs = np.abs(tmpDiff)
        # tmpDiff = tmpDiff/len(tmpPilotLocation)
        tmpDiff = tmpDiff/tmpAbs

        tmpOut = [each / tmpDiff  for each in inQam ]
        self.pilotPIdx = (self.pilotPIdx  + 1) % 127
        
        return tmpOut

    def htVhtCrc8Check(self,inBits): # in 2 * 24 bits. check B10-B19 in ht-sig2 and vht-sig-A2
        rxCrc8 = inBits[34:42]
        calCrc8 = p8h.genBitBitCrc8(inBits[0:34])
        if rxCrc8 == calCrc8:
            return True
        else: return False

    def proRxSigSym(self,inSym,mode): # inSym 80 samples  #L-SIG, HT-SIG, VHT-SIGA
        tmp = []
        # compensate cfo 
        print("self.symProcIdx in procRxSig, self.nSymProcIdx and self.pilotPIdx ",self.nSymProcIdx,self.pilotPIdx )
        
        for i in range(len(inSym)):
            tmp.append(inSym[i]*complex(np.cos(i*self.cfoRad+(self.nSymProcIdx*80)), np.sin(i*self.cfoRad+ (self.nSymProcIdx*80)))) 
        # remove CP 
        tmpTimeSym = tmp[16:80]

        #convert to freq domain and process symbol
        tmpSigFreq = p8h.procFftDemod(tmpTimeSym)
        tmpSigFreq = self.procCompChan(tmpSigFreq)
        tmpSigFreq = p8h.procRmDcNonDataSc(tmpSigFreq, p8h.F.L)
        tmpSigFreq = self.procPilotTrack(tmpSigFreq,p8h.C_PILOT_L)
        tmpSigLlr = p8h.procRemovePilots(tmpSigFreq)
        
        
        
        if mode == p8h.M.BPSK:
            tmpSigCoded = list(np.real(p8h.procDeinterleaveSigL(tmpSigLlr)))
        elif mode == p8h.M.QBPSK:
            tmpSigCoded = list(np.imag(p8h.procDeinterleaveSigL(tmpSigLlr)))

        self.nSymProcIdx += 1
        return tmpSigCoded
    
    def __procRxLegacySigSym(self): #inSym in time domian, 80 smaples
        tmpInSym = self.timeInSig[self.legacyStfIndex+320:self.legacyStfIndex+320+80] 
        
        self.LSigBits =  p8h.procViterbiDecoder(self.proRxSigSym(tmpInSym,p8h.M.BPSK),24, p8h.CR.CR12)
        print( self.LSigBits)


        tmpRecvHeaderBits = self.LSigBits[0:17]
        tmpParityBit = self.LSigBits[17]
        if sum(tmpRecvHeaderBits)%2 == tmpParityBit:
            print("Pass L-SIG legacy check")
            if self.ifdb:print("DEBUG: legacy sig bits lenght:",len(self.LSigBits))
            if self.ifdb:print(self.LSigBits)

            self.mcs = p8h.C_LEGACY_RATE_BIT.index(self.LSigBits[0:4])
            self.mdpuLen = self.bi2Deci(self.LSigBits[5:17])
            if self.ifdb:print("DEBUG: legacy length:",self.mdpuLen)
            if self.ifdb:print("DEBUG: legacy mcs:",self.mcs)
        else:    
            print("ERROR: legacy sig Parity Bits error")
            self.mcs = None
            return
        
       

    def procCompChan(self,fftSymIn):#in 64 samples
        # compensate channel
        tmpOut = []
        for i in range(0,64):
            if(self.chanInfoLtf[i] == 0):
                tmpOut.append(0)
            else:
                tmpOut.append(fftSymIn[i]/self.chanInfoLtf[i])
        return list(tmpOut)

    def bi2Deci(self,inArray):
        tmpOut = 0
        for i in range(len(inArray)):
            if inArray[i] == 1:
                tmpOut += 1<<i
        return tmpOut 


class dephy80211sumimo():
    def __init__(self,addr1,addr2,addr3="",addr4=""):
        self.addresses = [addr1,addr2,addr3,addr4] 
        self.stream = [] 
        self.d = []
        
        for addr in self.addresses:
            if addr:  
                tmpSisoPhy = dephy80211siso(addr,ifaddNoise = True ,ifDebug = True )
                self.d.append(tmpSisoPhy)#object for each stream, demod is created from siso
                self.stream.append(tmpSisoPhy.rxSisoStepsList()) 
        self.nSymIdx = 0
        self.nSS = len(self.stream)
        self.nLtf = self.nSS
        self.nSts =  self.nSS #shouldn't be the same. 
        self.mimonSym = self.d[0].demod.nSym #asuuming symbol number is equal in each stream 

        #channel 
        self.NLchanInv = []

        #data
        self.ssSymChanComp = [[] for i in range(self.nSS)]
        self.ssDataQam = [[] for i in range(self.nSS)]



    def rxMimoStepsList(self):
        print("in rxMimoStepsList",self.d[0].demod.nSym)
        self.__rxMimoNLchanEstimate()
        self.__rxMimoNLchanComp()
        self.__rxMimoProcData()
    def __rxMimoProcData(self):
          
        # tmpSigFreq = p8h.procFftDemod(tmpSymTime)
        # tmpSigFreq = self.procCompChan(tmpSigFreq) #64


        # tmpSigFreq = p8h.procRmDcNonDataSc(tmpSigFreq, self.demod.phyFormat)
        # tmpSigFreq = self.procPilotTrack(tmpSigFreq,tmpPilot)
        # print("len procPilotTrack", len(tmpSigFreq))
        # self.rxDataQam += p8h.procRemovePilots(tmpSigFreq) 
        # print("len  self.rxDataQam", len(self.rxDataQam))
        # self.nSymProcIdx +=1 
        # if(not(self.phyFormat == p8h.F.L)):
        #     tmpPilot = tmpPilot[1:] + [tmpPilot[0]] #?
        # # myConstellationPlot( self.rxDataQam)
        tmpPilot = []
        tmpPilotPIdx = 3 #data symbol pilot index starting at 3
        print("self.d[ssIter].demod",self.d[0].pilotPIdx)

        for ssIter in range(1):#self.nSS
            tmpPilot.append (p8h.C_PILOT_HT[p8h.BW.BW20.value][self.nSS - 1][ssIter])

            self.d[ssIter].pilotPIdx = 3 
            for symIter in range(1):#self.mimonSym
                tmpSym =  self.ssSymChanComp[ssIter][symIter*64:(symIter+1)*64]
                print("len",len(tmpSym))
                tmpSym = p8h.procRmDcNonDataSc(tmpSym, self.d[ssIter].demod.phyFormat)
                print("len",len(tmpSym))
                myConstellationPlot(tmpSym)
                tmpSym = self.d[ssIter].procPilotTrack(tmpSym,tmpPilot[ssIter])
                myConstellationPlot(tmpSym)
                if(not(self.d[ssIter].demod.phyFormat == p8h.F.L)):
                    tmpPilot[ssIter] = tmpPilot[ssIter][1:] + [tmpPilot[ssIter][0]] 

            print("----")
    def __rxMimoNLchanComp(self):
      
        # dataSym = self.d[0].demod.nSym
        for symIter in range(1):#self.mimonSym
            

            tmpSS= [] #one symbol
            #fft first
            for ssItr in range(0, self.nSS): #len of tmpSS = nSS, taking first symbol
                tmp = self.stream[ssItr][(symIter+self.nSymIdx)*64:(symIter+1+self.nSymIdx)*64]
                tmpFreq = p8h.procFftDemod(tmp)
                myConstellationPlot(tmpFreq)
                tmpSS.append(tmpFreq)
            tmp64sc = [] #([[ss1d1],[ss2d1]], [[ss1d2],[ss2d2]])
                    # print(np.matmul(rxltfSc,htChanInv))#to recover TX 

            for scIter in range(64):
                tmp1sc = []
                for ssItr in range(0, self.nSS):
                    tmp1sc.append([tmpSS[ssItr][scIter]])
                tmp64sc.append(tmp1sc)
            
                tmp1compOut = np.matmul(self.NLchanInv[scIter], tmp1sc )
                    # print(np.matmul(rxltfSc,htChanInv))#to recover TX 

                #sidtribute comp out to stream 
                for ssItr in range(0, self.nSS):
                     self.ssSymChanComp[ssItr].append(tmp1compOut[ssItr])
   
            print("  self.ssSymChanComp",len(  self.ssSymChanComp),len(  self.ssSymChanComp[0]),  self.ssSymChanComp)
    def __rxMimoNLchanEstimate(self):
        self.nSymIdx = 1 #skip HT-STF 
        tmpChanSyms = [] 
        htLtfTxP = []
        htChanInv = []
        for ssItr in range(0, self.nSS): #nLTF = tx number   #tmpChanSyms has nSTS(rx antenna) number    
            rxltfSym = [] #from received streams
            tmpLtfP = [] 
            rxltfSc = []
            
            for ltfIter in range(0, self.nLtf):#[56*2,56*2]
                tmpLtfP.append(p8h.C_P_LTF_HT_4[ssItr][ltfIter])
                rxltfSym += list(p8h.procFftDemod(self.stream[ssItr][(self.nSymIdx+ltfIter)*64:(self.nSymIdx+ltfIter)*64+64]))
                print(len(rxltfSym))

            htLtfTxP.append(tmpLtfP)
            tmpChanSyms.append(rxltfSym)
        for scIter in range(64): #to make tmpStsByLtf has 64 2by2 rxSym
            tmpStsByLtf = []
            for ssItr in range(0, self.nSS):
                tmpLtf = [] 
                for ltfIter in range(0, self.nLtf):
                    tmpLtf.append(self.stream[ssItr][64*ltfIter+scIter])
                tmpStsByLtf.append(tmpLtf)
            rxltfSc.append(np.array(tmpStsByLtf))

            tmpChan = np.matmul(np.linalg.inv(htLtfTxP),tmpStsByLtf)
            htChanInv.append(np.linalg.inv(tmpChan))
        self.NLchanInv  = htChanInv
        # print(np.matmul(rxltfSc,htChanInv))#to recover TX 

        self.nSymIdx+=self.nLtf # for ht, after ht-LTF is data symbol




# p8h.C_SCALENTF_LTF_VHT[self.m.bw.value]


        # myConstellationPlot(tmpChanSyms[0][0:56])
        # myConstellationPlot(tmpChanSyms[1][0:56])
        # myConstellationPlot(tmpChanSyms[1][0+56:56+56])
def myConstellationPlot(inSig):
        x = []
        y = [] 
        absArray = []
        for i in range(len(inSig)):
            x.append(np.real(inSig[i]))
            y.append(np.imag(inSig[i]))
            absArray.append(cmath.polar(inSig[i]))
        # print(absArray)
        plt.title("ConstellationPlot in self")
        # plt.xlim([-1.1, 1.1])
        # plt.ylim([-1.1, 1.1])
        # plt.xlim([-2, 2])
        # plt.ylim([-2, 2])
        plt.scatter(x,y)
        plt.show()

def myPlot(inSigComp):
    plt.plot(np.real(inSigComp))
    plt.plot(np.imag(inSigComp))
    plt.show()

def myStemPlot(inSigCplx):
    x = []
    y = [] 
    for i in range(len(inSigCplx)):
        x.append(np.real(inSigCplx[i]))
        y.append(np.imag(inSigCplx[i]))
    plt.title("ConstellationPlot in self")
    # plt.xlim([-1.1, 1.1])
    # plt.ylim([-1.1, 1.1])
    plt.stem(x,y)
    plt.show()

def compare(in1, in2):
    
    count = 0
    tmpOut = True
    for i in range(min(len(in1),len(in2))):
        if in1[i]!=in2[i]:
            # print("index",i)
            count+=1
            tmpOut = False
            
    
    return (count,tmpOut)


def procScramble(inBits, scrambler):
    tmpScrambleBits = []
    for i in range(0, len(inBits)):
        tmpFeedback = int(not not(scrambler & 64)) ^ int(not not(scrambler & 8))
        tmpScrambleBits.append(tmpFeedback ^ inBits[i])
        scrambler = ((scrambler << 1) & 0x7e) | tmpFeedback
    return tmpScrambleBits

def biArray2PackByte(inBinary):
    tmpOut = b''

    for i in range(int(len(inBinary)/8)):
        tmpBi = inBinary[i*8:i*8+8]

        tmpSum = 0
        for i in range(8):
            if tmpBi[i] == 1:
                tmpSum += 1<<i
        tmpOut += struct.pack('>B',tmpSum)

    return tmpOut


if __name__ == "__main__":



    pyToolPath = os.path.dirname(__file__)
    "siso rx "
    # addr = os.path.join(pyToolPath, "../tmp/sig80211GenMultipleSiso_1x1_0.bin")
    # # addr = os.path.join(pyToolPath, "../tmp/recordchop.bin")
    # phy = dephy80211siso(addr,ifaddNoise = False ,ifDebug = True )
    # phyBits = phy.rxSisoStepsList()

    # phyBiToBytePack = biArray2PackByte(phyBits) # to match org generated in TX 
    # print("hex compare with mac pkt ->", phyBiToBytePack == b'\x08\x01n\x00\xf4i\xd5\x80\x0f\xa0\x00\xc0\xca\xb1[\xe1\xf4i\xd5\x80\x0f\xa0\x00\xa9\xaa\xaa\x03\x00\x00\x00\x08\x00E\x00\x00:\xab\x02@\x00@\x11{\x96\n\n\x00\x06\n\n\x00\x01\x99\xd3"\xb9\x00&\x10\xec123456789012345678901234567890\xa3]\xee\xec')

    "mimo rx"
    # addr1 = os.path.join(pyToolPath,"../tmp/sig80211GenMultipleMimo_2x2_0.bin")
    # addr2 = os.path.join(pyToolPath,"../tmp/sig80211GenMultipleMimo_2x2_1.bin")
    addr1 = os.path.join(pyToolPath,"../tmp/ConvOut0.bin")
    addr2 = os.path.join(pyToolPath,"../tmp/ConvOut1.bin")
    phymimo = dephy80211sumimo(addr1,addr2)
    phymimo.rxMimoStepsList()

    

    # MAC       
    # upperlayer = mac.RxMac80211(phyBits, debug = True)
    # msg = upperlayer.rxMacStepsList()
    # print("udp payload:",msg)
    
    
    