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
class dephy80211():
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
        self.phyFormate = None

        #rx phy bits out 
        self.phyBits = []
        # debug
        self.ifdb = ifDebug

    def rxStepsList(self): 
        # self.__testNoCfo()
        self.__procRxLegacyStfTrigger()
        self.__procRxLegacyCfoEst()
        self.__procRxLegacyChanEst()
        self.__procRxLegacySigSym()
        self.__procRxPktFormat()
        if self.nSS == 1:
            if self.phyFormat == p8h.F.L:
                self.__procRxSisoData()
            elif self.phyFormat == p8h.F.VHT:
                self.__procRxChanUpdate
                self.__procRxParseVHTSigB()
                self.__procRxSisoData()
            elif self.phyFormat == p8h.F.HT:
                self.__procRxChanUpdate()
                self.__procRxSisoData()
        else:
            pass # more than one ht-LTF

        


        




        return self.rxDatabits
    

    def __procRxParseVHTSigB(self):
        #need to decode VHT-SIGB
        pass


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
                                    # after channel compensate, processing symbol
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
        print("channel updated")

    def __procRxSisoData(self):
        if self.ifdb:print("DEBUG: final format:",self.phyFormat)
        if self.ifdb:print("DEBUG: final mcs:",self.mcs)        
        self.demod = p8h.modulation(self.phyFormat, self.mcs, bw=p8h.BW.BW20, nSTS=self.nSS, shortGi=False)
        self.demod.procPktLenNonAggre(self.mdpuLen)

        print("self.nSymProcIdx",self.nSymProcIdx) #siso doesn't utilize HT-STF and HT LTF(s)
        if self.phyFormat == p8h.F.L:
            self.nSymProcIdx = 5
        elif  self.phyFormat == p8h.F.HT:
            self.nSymProcIdx = 9 
        elif self.phyFormat == p8h.F.VHT:
            self.nSymProcIdx = 10

        self.rxDataSym = self.timeInSig[self.legacyStfIndex+(self.nSymProcIdx*80):self.legacyStfIndex+(self.nSymProcIdx*80)+self.demod.nSym*80]
        print("data sample length",len(self.rxDataSym), self.demod.nSym)
        # myPlot(self.rxDataSym)
        # self.demod.nSym = 1
        tmpPilot = p8h.C_PILOT_L
        for nSymIter in range(self.demod.nSym):#self.demod.nSym
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
            print("len procPilotTrack", len(tmpSigFreq))
            self.rxDataQam += p8h.procRemovePilots(tmpSigFreq) 
            print("len  self.rxDataQam", len(self.rxDataQam))
            self.nSymProcIdx +=1 
            if(not(self.phyFormat == p8h.F.L)):
                tmpPilot = tmpPilot[1:] + [tmpPilot[0]] #?
            # myConstellationPlot( self.rxDataQam)



        if self.ifdb:print("DEBUG: self.rxDataQam len:",len(self.rxDataQam))
        # if self.ifdb:print("DEBUG: self.rxDataQam", self.rxDataQam)
        self.rxInterleaveBits = self.procSymQamToLlr()
        print("compare interleave",compare([1 if np.real(each)>0 else 0  for each in self.rxInterleaveBits ],txbits.htMcs0interBits))
        # if self.ifdb:print("DEBUG: self.rxInterleaveBits len:",len(self.rxInterleaveBits))
        # if self.ifdb:print("DEBUG: self.rxInterleaveBits: \n", self.rxInterleaveBits)

        self.rxCodedBits = self.procDeinterleave()
        if self.ifdb:print("DEBUG: self.rxCodedBits len:",len(self.rxCodedBits))
        print("compare coded",compare( [1 if np.real(each)>0 else 0  for each in self.rxCodedBits ],txbits.htMcs0CodedBits))

        
        # if self.ifdb:print("DEBUG: self.rxCodedBits: \n", self.rxCodedBits)
        
        self.rxScrambleBits = p8h.procViterbiDecoder(self.rxCodedBits,self.demod.nDBPS*self.demod.nSym,self.demod.cr)
        if self.ifdb:print("DEBUG: self.rxScrambleBits len:",len(self.rxScrambleBits))
        # if self.ifdb:print( "DEBUG: self.rxScrambleBits",self.rxScrambleBits)

        self.rxDatabits = self.PrcoDeScrambleBits()
        if self.ifdb:print("DEBUG: self.rxDatabits len:", len(self.rxDatabits) )
        # if self.ifdb:print("DEBUG: self.rxDatabits \n", self.rxDatabits )
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

        return tmpDeScrambledBits[16:self.mdpuLen * 8 +16]#removing 16 serving bits and 6 tail bits  

    def procDeinterleave(self):
        if self.phyFormat == p8h.F.L:
            tmpDeinterleaveBits = [0] * self.demod.nSym * self.demod.nCBPS
            print("in procDeinterleave",self.demod.phyFormat, self.demod.cr, self.demod.nSym * self.demod.nCBPS, self.demod.nCBPS)
            s = int(max(1, self.demod.nBPSCS/2))
            for symIter in range(self.demod.nSym):
                for j in range(self.demod.nCBPS):
                    i = int( int(s * np.floor(j/s)) +  int((j+np.floor(16*j/self.demod.nCBPS))%s))  
                    k = int(16 * i - (self.demod.nCBPS - 1) * int(np.floor(16 * i /self.demod.nCBPS)))
                    tmpDeinterleaveBits[symIter*self.demod.nCBPS+k] =  self.rxInterleaveBits[symIter*self.demod.nCBPS+j]
            
        elif self.phyFormat== p8h.F.HT:
            print("in procDeinterleave ht ")
            tmpDeinterleaveBits =  [0] * self.demod.nSym * self.demod.nCBPSS
            s = int(max(1, self.demod.nBPSCS/2))
            for symIter in range(self.demod.nSym):#self.demod.nSym
                if self.nSS == 1:
                    for r in range(self.demod.nCBPSS):
                        j = r
                        i = s * int(np.floor(j/s)) + (j+ int(np.floor(self.demod.nIntlevCol*j/self.demod.nCBPSS)))%s
                        k = self.demod.nIntlevCol * i - (self.demod.nCBPSS - 1 ) * int( np.floor(i / self.demod.nIntlevRow))
                        tmpDeinterleaveBits[symIter*self.demod.nCBPSS+k] = self.rxInterleaveBits[symIter*self.demod.nCBPSS+r]
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
    
    def __procRxPktFormat(self):#here
        tmpSig1 = self.timeInSig[self.legacyStfIndex+400:self.legacyStfIndex+400+80]
        tmpSig2 = self.timeInSig[self.legacyStfIndex+480:self.legacyStfIndex+480+80]
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
                    print("htCrc8Check pass")
                    print("len HT", self.bi2Deci(tmpHtSigBits[8:24]))
                    self.mdpuLen = self.bi2Deci(tmpHtSigBits[8:24])
                    self.phyFormat = p8h.F.HT 
                    self.mcs = self.bi2Deci(tmpHtSigBits[0:7])

                if (self.htVhtCrc8Check(tmpVhtSigBits)):
                    if self.ifdb:print("DEBUG: VHT-SIGA  bits:",tmpHtSigBits)
                    print("vhtCrc8Check pass")
                    self.pilotPIdx = 4
                    self.phyFormat = p8h.F.VHT
                    self.mcs = self.bi2Deci(tmpVhtSigBits[28:32])

                if (not (self.htVhtCrc8Check(tmpHtSigBits) or self.htVhtCrc8Check(tmpVhtSigBits))):
                    print("state machin to L mcs 0 ")
                    self.pilotPIdx = 1 
                    self.phyFormat = p8h.F.L
            else: #demodulate legacy mcs > 0 
                #??how to chexk bandwidth parameter 
                self.phyFormat = p8h.F.L
                self.pilotPIdx = 1 
                
            self.nSS = int(np.floor(self.mcs / 8)) + 1
        else:
            return "ERROR: legacy Sig parity check fail"
            
    def procPilotTrack(self,inQam,IN_PILOT):
        tmpOut = None
        tmpPilotLocation = None
        tmpStandardPilot = [int (each * p8h.C_PILOT_PS[self.pilotPIdx]) for each in IN_PILOT ]
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

    def proRxSigSym(self,inSym,mode): # inSym 80 samples 
        outBits = None
        tmp = []
        # compensate cfo 
        print("self.symProcIdx in procRxSig",self.nSymProcIdx)
        
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
        # myConstellationPlot(tmpSigFreq)
        
        
        
        if mode == p8h.M.BPSK:
            tmpSigCoded = list(np.real(p8h.procDeinterleaveSigL(tmpSigLlr)))
        elif mode == p8h.M.QBPSK:
            tmpSigCoded = list(np.imag(p8h.procDeinterleaveSigL(tmpSigLlr)))

        self.nSymProcIdx += 1
        return tmpSigCoded
    
    def __procRxLegacySigSym(self): #inSym in time domian, 80 smaples
        tmpInSym = self.timeInSig[self.legacyStfIndex+320:self.legacyStfIndex+320+80] 
        
        self.LSigBits =  p8h.procViterbiDecoder(self.proRxSigSym(tmpInSym,p8h.M.BPSK),24, p8h.CR.CR12)



        tmpRecvHeaderBits = self.LSigBits[0:17]
        tmpParityBit = self.LSigBits[17]
        if sum(tmpRecvHeaderBits)%2 == tmpParityBit:
            print("Pass L-SIG legacy check")
            if self.ifdb:print("DEBUG: legacy sig bits lenght:",len(self.LSigBits))
            if self.ifdb:print(self.LSigBits)

            self.mcs = p8h.C_LEGACY_RATE_BIT.index(self.LSigBits[0:4])
            self.mdpuLen = self.bi2Deci(self.LSigBits[5:17])
            if self.ifdb:print("DEBUG: legacy sig mpdu length:",self.mdpuLen)
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
        plt.xlim([-2, 2])
        plt.ylim([-2, 2])
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
    addr = os.path.join(pyToolPath, "../tmp/sig80211GenMultipleSiso_1x1_0.bin")
    phy = dephy80211(addr,ifaddNoise = False ,ifDebug = True )
    phyBits = phy.rxStepsList()

    # phyBiToBytePack = biArray2PackByte(phyBits) # to match org generated in TX 


    # MAC       
    
    # phyBits = [1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1]
    # upperlayer = mac.RxMac80211(phyBits, debug = True)
    # msg = upperlayer.rxMacStepsList()
    # print("udp payload:",msg)
    
    
    
    #tx ht sig
    #[0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    #[0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    
    #VHT SIG A mcs1
    #[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    

    # input = [-1.0444180550565167, -0.9947426356225831, -1.0401742178810855, -1.031713245740612, 1.0370877177155429, 0.9913488535371258, -0.9888611001366912, 1.0350719123748728, 1.0241458393922958, 1.013705415285691, 1.0455470602091317, 1.0437183787661155, -1.023272836151735, -1.0325727413044568, 1.015976243553799, -0.9996854183591589, 0.9952186351218241, 1.0077086309153938, -0.977323838476975, -1.0281748240680866, -1.0159761901837367, -1.0024165611893632, -0.9611991904296747, -0.9729657812260926, -0.9764090775993102, -0.9570488677666111, -0.9564583109545355, -1.015800732350474, -0.9667693199912082, -1.0350181092291935, -1.0408546761550646, -1.0039115477276974, -0.9985348557220102, -0.9838687472656028, -1.0349162084956103, -1.0243422431655052, -0.9618928728978773, -0.9978311163549579, -0.980464476975575, -0.9696386500576356, -0.9797858905983382, -0.9911452419586753, -1.0077085980739628, -1.0196368231112465, -0.9696116108768597, -0.9796382879214189, 0.9602369107790859, 0.9604122669883367,-1.0451372819658478, 0.9367985018516418, 0.9804009406325066, 0.9787893806065077, 0.9791310425454565, 0.9349056125167476, -1.041086094375738, -1.0185353898205634, -0.9871579046266233, 1.0057449478757279, 1.0579690583980519, -0.9824829813661166, 1.08100180514853, 1.0815394415078197, 1.0743864337353985, 1.0152352820884982, -1.051797883094749, -1.0685512247040998, -1.018355715572064, 0.9786841225320102, 0.9796761569749066, -1.0603788260359572, 0.9623195330104091, 0.9941083499298758, -0.9982835380364837, -0.9586429111149045, 1.0136861829357897, -0.9813481975427704, -0.9218506267778753, -0.9787987407862584, 1.0657250863459635, 1.0140265302206024, 0.9414717425930793, 0.9263920305955843, 0.979113482717666, -1.0799858708751677, -1.0179398543557732, 0.9837949183728565, 0.9903142001665646, 0.9960693125930508, -0.9962709192704778, -0.9904095027225215, -0.9837910955713197, -0.9804009348933734, -0.9216969338563331, -0.9253324765844163, -0.9393278459040605, -0.9558574020548589] 
    # outBits = p8h.procViterbiDecoder(input, 48, p8h.CR.CR12)
    # print("output",outBits )