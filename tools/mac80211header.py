"""
    GNU Radio IEEE 802.11a/g/n/ac 2x2
    Python tools
    Copyright (C) June 1, 2022  Zelin Yun

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import socket
import mac80211
import phy80211header as p8h
from matplotlib import pyplot as plt
import numpy as np
import time
import struct
from enum import Enum

# number of bit for each angle to quantize for V of 2 rows, 3 rows and 4 rows
C_VHT_CB0_R2_ANGLE_BIT_N = [7, 5]
C_VHT_CB0_R3_ANGLE_BIT_N = [7, 7, 5, 5, 7, 5]
C_VHT_CB0_R4_ANGLE_BIT_N = [7, 7, 7, 5, 5, 5, 7, 7, 5, 5, 7, 5]
C_VHT_CB1_R2_ANGLE_BIT_N = [9, 7]
C_VHT_CB1_R3_ANGLE_BIT_N = [9, 9, 7, 7, 9, 7]
C_VHT_CB1_R4_ANGLE_BIT_N = [9, 9, 9, 7, 7, 7, 9, 9, 7, 7, 9, 7]

class FC_TPYE(Enum):
    MGMT = 0
    CTRL = 1
    DATA = 2
    EXT = 3

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

class FC_SUBTPYE_MGMT(Enum):
    ASSOREQ = 0
    ASSORES = 1
    REASSOREQ = 2
    REASSORES = 3
    PROBEREQ = 4
    PROBERES = 5
    TIMINGAD = 6
    RESERVED7 = 7
    BEACON = 8
    ATIM = 9
    DISASSO = 10
    AUTH = 11
    DEAUTH = 12
    ACT = 13
    ACTNOACK = 14
    RESERVED15 = 15

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

class FC_SUBTPYE_CTRL(Enum):
    RESERVED0 = 0
    RESERVED1 = 1
    RESERVED2 = 2
    RESERVED3 = 3
    BFREPOPOLL = 4
    VHTNDPANNO = 5
    FRAMEEXT = 6
    WRAPPER = 7
    BLOCKACKREQ = 8
    BLOCKACK = 9
    PSPOLL = 10
    RTS = 11
    CTS = 12
    ACK = 13
    CFEND = 14
    CFENDCFACK = 15

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

class FC_SUBTPYE_DATA(Enum):
    DATA = 0
    DATACFACK = 1
    DATACFPOLL = 2
    DATACFACKCFPOLL = 3
    NULL = 4
    CFACK = 5
    CFPOLL = 6
    CFACKCFPOLL = 7
    QOSDATA = 8
    QOSDATACFACK = 9
    QOSDATACFPOLL = 10
    QOSDATACFACKCFPOLL = 11
    QOSNULL = 12
    RESERVED13 = 13
    QOSCFPOLL = 14
    QOSCFACKCFPOLL = 15

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

class FC_SUBTPYE_EXT(Enum):
    DMGBEACON = 0
    RESERVED1 = 1
    RESERVED2 = 2
    RESERVED3 = 3
    RESERVED4 = 4
    RESERVED5 = 5
    RESERVED6 = 6
    RESERVED7 = 7
    RESERVED8 = 8
    RESERVED9 = 9
    RESERVED10 = 10
    RESERVED11 = 11
    RESERVED12 = 12
    RESERVED13 = 13
    RESERVED14 = 14
    RESERVED15 = 15

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

class MGMT_ELE(Enum):
    SSID = 0
    SUPOTRATE = 1
    DSSSPARAM = 3
    TIM = 5
    COUNTRY = 7
    BSSLOAD = 11
    HTCAP = 45
    RSN = 48
    HTOPS = 61
    ANTENNA = 64
    RMENABLED = 70
    EXTCAP = 127
    VHTCAP = 191
    VHTOPS = 192
    TXPOWER = 195
    VENDOR = 221

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

# some define
C_FC_SUBTYPE_MGMT_STR = ["Association Request", "Association Response", "Reassociation Request", "Reassociation Response", "Probe Request", "Probe Response", "Timing Advertisement", "Reserved", "Beacon", "ATIM", "Disassociation", "Authentication", "Deauthentication", "Action", "Action No Ack", "Reserved"]
C_FC_SUBTYPE_CTRL_STR = ["Reserved", "Reserved", "Reserved", "Reserved", "Beamforming Report Poll", "VHT NDP Announcement", "Control Frame Extension", "Control Wrapper", "Block Ack Request (BlockAckReq)", "Block Ack (BlockAck)", "PS-Poll", "RTS", "CTS", "Ack", "CF-End", "CF-End +CF-Ack"]
C_FC_SUBTYPE_DATA_STR = ["Data", "Data +CF-Ack", "Data +CF-Poll", "Data +CF-Ack +CF-Poll", "Null (no data)", "CF-Ack (no data)", "CF-Poll (no data)", "CF-Ack +CF-Poll (no data)", "QoS Data", "QoS Data +CF-Ack", "QoS Data +CF-Poll", "QoS Data +CF-Ack +CF-Poll", "QoS Null (no data)", "Reserved", "QoS CF-Poll (no data)", "QoS CF-Ack +CF-Poll (no data)"]
C_FC_SUBTYPE_EXT_STR = ["DMG Beacon", "Reserved", "Reserved", "Reserved", "Reserved","Reserved", "Reserved", "Reserved", "Reserved","Reserved", "Reserved", "Reserved", "Reserved", "Reserved", "Reserved", "Reserved"]

class frameControl:
    def __init__(self, fc):
        self.fcValue = fc
        self.protocalVer = fc & 3
        self.type = FC_TPYE((fc >> 2) & 3)
        if(self.type == FC_TPYE.MGMT):
            self.subType = FC_SUBTPYE_MGMT((fc >> 4) & 15)
        elif(self.type == FC_TPYE.CTRL):
            self.subType = FC_SUBTPYE_CTRL((fc >> 4) & 15)
        elif(self.type == FC_TPYE.DATA):
            self.subType = FC_SUBTPYE_DATA((fc >> 4) & 15)
        elif(self.type == FC_TPYE.EXT):
            self.subType = FC_SUBTPYE_EXT((fc >> 4) & 15)
        else:
            self.subType = str((fc >> 4) & 15)
        self.toDs = (fc >> 8) & 1
        self.fromDs = (fc >> 9) & 1
        self.moreFrag = (fc >> 10) & 1
        self.retry = (fc >> 11) & 1
        self.powerMgmt = (fc >> 12) & 1
        self.moreData = (fc >> 13) & 1
        self.protectFrame = (fc >> 14) & 1
        self.htcOrder = (fc >> 15) & 1      # in non qos data set 1 for order, in qos data set 1 for htc
    
    def printInfo(self):
        print("cloud mac80211header, value %s, FC Info protocol:%d, type:%s, sub type:%s, to DS:%d, from DS:%d, more frag:%d, retry:%d" % (hex(self.fcValue), self.protocalVer, self.type, self.subType, self.toDs, self.fromDs, self.moreFrag, self.retry))

def procVhtPhiQuanti(phi, bphi):
    tmpK = list(range(0, 2**bphi))
    tmpShift = np.pi / (2**(bphi))
    tmpPhi = [each * np.pi / (2**(bphi-1)) + tmpShift for each in tmpK]
    # print(tmpPhi)
    return min(range(len(tmpPhi)), key=lambda i: abs(tmpPhi[i]-phi))

def procVhtPsiQuanti(psi, bpsi):
    tmpK = list(range(0, 2**bpsi))
    tmpShift = np.pi / (2**(bpsi+2))
    tmpPsi = [each * np.pi / (2**(bpsi+1)) + tmpShift for each in tmpK]
    # print(tmpPsi)
    return min(range(len(tmpPsi)), key=lambda i: abs(tmpPsi[i]-psi))

def procVhtVCompress(v, codeBookInfo = 0, ifDebug = False):
    resValue = []       # value and type
    resType = []
    # Phi for Di, Psi for Gli, feedback type is MU-MIMO for VHT PPDU
    if(isinstance(v, np.ndarray) and isinstance(ifDebug, bool) and isinstance(codeBookInfo, int)):
        [m, n] = v.shape    # m is row, n is col
        if(ifDebug):
            print("V, get V %d row %d col" % (m, n))
            print(v)
        if(m > 0 and n > 0 and m >= n):
            if(codeBookInfo):
                nBitPhi = 9
                nBitPsi = 7
            else:
                nBitPhi = 7
                nBitPsi = 5
            if(ifDebug):
                print("quantization codebook %d, Phi %d bits, Psi %d bits" % (codeBookInfo, nBitPhi, nBitPsi))
            dt = np.zeros((n, n), dtype=complex)
            for j in range(0, n):
                tmpTheta = np.arctan2(np.imag(v[m-1][j]), np.real(v[m-1][j]))
                tmpValue = np.exp(tmpTheta * 1.0j)
                dt[j][j] = tmpValue
            if(ifDebug):
                print("Dt, get D tilde")
                print(dt)
            dtH = dt.conjugate().T
            if(ifDebug):
                print("DtH, get D tilde hermitian")
                print(dtH)
            vdtH = np.matmul(v, dtH)
            if(ifDebug):
                vdtHRes = np.matmul(v, dtH)
            if(ifDebug):
                print("VDtH, get V dot D tilde hermitian, which is also the matrix with all real number on the last row, to be decomposed with givens rotation")
                print(vdtH)
            for j in range(0, n):
                vdtH[m-1][j] = np.real(vdtH[m-1][j])
            if(ifDebug):
                print("VDtH, get V dot D tilde hermitian, remove residual imag for last row")
                print(vdtH)
            
            glidiHvdtH_name = "VDtH"
            if(ifDebug):
                vtRes = np.identity(m)
            # each loop we compute Di and following Gli(s), and keep updating the vdtH to be GliDiHVDtH
            for grIter in range(0, min(m-1, n)):      # gr for givens rotation
                i = grIter + 1
                if(ifDebug):
                    print("givens rotation loop round %d" % (i))
                di = np.zeros((m, m), dtype=complex)
                for j in range(0, i-1):
                    di[j][j] = 1
                diPhi = []
                for j in range(i, m):
                    tmpTheta = np.arctan2(np.imag(vdtH[j-1][i-1]), np.real(vdtH[j-1][i-1]))
                    diPhi.append(tmpTheta)
                    if(ifDebug):
                        print("D%d Phi:%f" % (i, tmpTheta))
                    tmpValue = np.exp(tmpTheta * 1.0j)
                    di[j-1][j-1] = tmpValue
                if(len(diPhi)):
                    diPhi = np.unwrap(diPhi)
                    if(diPhi[0] < 0):
                        diPhi += np.pi * 2
                    for each in diPhi:
                        resValue.append(procVhtPhiQuanti(each, nBitPhi))
                        resType.append(0)
                    if(ifDebug):
                        print("D%d Unwrapped Phi:" % (i), diPhi)
                        print("D%d Unwrapped Phi Quantized:" % (i), [procVhtPhiQuanti(each, nBitPhi) for each in diPhi])
                di[m-1][m-1] = 1
                if(ifDebug):
                    print("D%d, get Di here" % (i))
                    print(di)
                if(ifDebug):                
                    vtRes = np.matmul(vtRes, di) # compute the final Vt to compare
                vdtH = np.matmul(di.conjugate().T, vdtH)    # now vdtH is from VDtH to DiHVDtH
                if(ifDebug):
                    glidiHvdtH_name = "D%dH" % (i) + glidiHvdtH_name
                    print(glidiHvdtH_name + ", get D%d hermitian dot V dot D tilde hermitian" % (i))
                    print(vdtH)
                # remove residual imag of the column
                for l in range(i, m):
                    vdtH[l-1][i-1] = np.real(vdtH[l-1][i-1])
                if(ifDebug):
                    print(glidiHvdtH_name + ", the %dth column should be all real" % (i))
                    print(vdtH)
                for l in range(i+1, m+1):
                    if(ifDebug):
                        print("compute givens rotation G%d%d" % (l, i))
                    gli = np.zeros((m, m), dtype=complex)
                    x1 = np.real(vdtH[i-1][i-1])
                    x2 = np.real(vdtH[l-1][i-1])
                    y = np.sqrt(x1*x1 + x2*x2)
                    gliPsi = np.arccos(x1 / y)
                    resValue.append(procVhtPsiQuanti(gliPsi, nBitPsi))
                    resType.append(1)
                    if(ifDebug):
                        print("x1:%f, x2:%f, y:%f, GliPsi:%f, Quantized GliPsi:%d" % (x1, x2, y, gliPsi, procVhtPsiQuanti(gliPsi, nBitPsi)))
                    if(gliPsi < 0 or gliPsi > np.pi/2):
                        print("Gli psi value error!!!!!!!!!!!!!!!!!!!!")
                    for j in range(0, m):
                        gli[j][j] = 1
                    gli[i-1][i-1] = np.cos(gliPsi)
                    gli[l-1][i-1] = -np.sin(gliPsi)
                    gli[i-1][l-1] = np.sin(gliPsi)
                    gli[l-1][l-1] = np.cos(gliPsi)
                    if(ifDebug):
                        print("G%d%d" % (l, i))
                        print(gli)
                    vdtH = np.matmul(gli, vdtH)
                    gliT = gli.T
                    if(ifDebug):
                        vtRes = np.matmul(vtRes, gliT) # compute the final Vt to compare
                    if(ifDebug):
                        glidiHvdtH_name = "G%d%d" % (l, i) + glidiHvdtH_name
                        print(glidiHvdtH_name + ", now the %d%d location shoule be zero" % (l, i))
                        print(vdtH)
                    vdtH[l-1][i-1] = 0
                    if(ifDebug):
                        print(glidiHvdtH_name + ", remove %d%d residual error" % (l, i))
                        print(vdtH)
            if(ifDebug):
                vIt = np.zeros((m, n), dtype=complex) # I tilde for V tilde
                for j in range(0, min(m, n)):
                    vIt[j][j] = 1
                vtRes = np.matmul(vtRes, vIt)
                print("compare the VDtH and decompsed results")
                print(vdtHRes)
                print(vtRes)
                print(resValue)
                print(resType)
    return resValue, resType

def mgmtElementParser(inbytes):
    if(isinstance(inbytes, (bytes, bytearray)) and len(inbytes) > 0):
        elementIter = 0
        tmpMgmtElements = []
        inPutLen = len(inbytes)
        while((elementIter+2) < inPutLen):  # one byte type, one byte len
            # print("Element at %d, ID %d, Len %d" % (elementIter, inbytes[elementIter], inbytes[elementIter+1]))
            if(MGMT_ELE.has_value(inbytes[elementIter])):
                tmpElement = MGMT_ELE(inbytes[elementIter])
                elementIter += 1
                tmpLen = inbytes[elementIter]
                elementIter += 1
                if((elementIter+tmpLen) < inPutLen):
                    if(tmpElement == MGMT_ELE.SSID):
                        print(inbytes[elementIter:elementIter+tmpLen].hex())
                        tmpStr = inbytes[elementIter:elementIter+tmpLen].decode("utf-8")
                        tmpMgmtElements.append("SSID: " + tmpStr)
                    elif(tmpElement == MGMT_ELE.SUPOTRATE):
                        tmpStr = ""
                        for i in range(0, tmpLen):
                            if(inbytes[elementIter+i] >= 0x80):
                                tmpStr += str((inbytes[elementIter+i]-0x80)*500/1000)
                                tmpStr += "Mbps(Basic) "
                            else:
                                tmpStr += str(inbytes[elementIter+i]*500/1000)
                                tmpStr += "Mbps "
                        tmpMgmtElements.append("Suppoprted Rates: " + tmpStr)
                    elif(tmpElement == MGMT_ELE.DSSSPARAM):
                        tmpStr = str(inbytes[elementIter])
                        tmpMgmtElements.append("DS Channel: " + tmpStr)
                    elif(tmpElement == MGMT_ELE.TIM):
                        tmpMgmtElements.append("TIM: To be added")
                    elif(tmpElement == MGMT_ELE.COUNTRY):
                        tmpStr = inbytes[elementIter:elementIter+3].decode("utf-8")
                        tmpMgmtElements.append("Country: " + tmpStr)
                    elif(tmpElement == MGMT_ELE.BSSLOAD):
                        tmpStaCount = struct.unpack('<H', inbytes[elementIter:elementIter+2])[0]
                        tmpChanUtil = int(inbytes[elementIter+2])
                        tmpAvailAdmCap = struct.unpack('<H', inbytes[elementIter+3:elementIter+5])[0]
                        tmpMgmtElements.append("BSS Load, Station Count: " + str(tmpStaCount) + ", Channel Utilization: " + str(tmpChanUtil) + ", Available Admission Capacity: " + str(tmpAvailAdmCap))
                    elif(tmpElement == MGMT_ELE.RSN):
                        tmpMgmtElements.append("RSN: To be added")
                    elif(tmpElement == MGMT_ELE.HTCAP):
                        tmpHtCapInfo = struct.unpack('<H', inbytes[elementIter:elementIter+2])[0]
                        tmpMgmtElements.append("HT CAP Info, LDPC %d, Chan Width %d, GF %d, SGI20 %d, SGI40 %d, TXSTBC %d" % (
                            (tmpHtCapInfo >> 0) & 1,
                            (tmpHtCapInfo >> 1) & 1,
                            (tmpHtCapInfo >> 4) & 1,
                            (tmpHtCapInfo >> 5) & 1,
                            (tmpHtCapInfo >> 6) & 1,
                            (tmpHtCapInfo >> 7) & 1
                        ))
                        tmpMcsBits = []
                        for i in range(0, 10):
                            for j in range(0, 8):
                                tmpMcsBits.append((inbytes[elementIter+3+i] >> j) & 1)
                        tmpMcsStr = "HT CAP MCS 0-31: "
                        for i in range(0, 32):
                            tmpMcsStr += str(tmpMcsBits[i])
                        tmpMgmtElements.append(tmpMcsStr)
                        tmpMcsStr = "HT CAP MCS 32-76: "
                        for i in range(32, 77):
                            tmpMcsStr += str(tmpMcsBits[i])
                        tmpMgmtElements.append(tmpMcsStr)
                    elif(tmpElement == MGMT_ELE.HTOPS):
                        tmpMgmtElements.append("HT Operation: To be added")
                    elif(tmpElement == MGMT_ELE.ANTENNA):
                        tmpMgmtElements.append("Antenna: %d" % inbytes[elementIter])
                    elif(tmpElement == MGMT_ELE.RMENABLED):
                        tmpMgmtElements.append("RM Enable Cap: To be added")
                    elif(tmpElement == MGMT_ELE.EXTCAP):
                        tmpMgmtElements.append("Ext Cap: To be added")
                    elif(tmpElement == MGMT_ELE.VHTCAP):
                        tmpVhtCapInfo = struct.unpack('<I', inbytes[elementIter:elementIter+4])[0]
                        tmpMgmtElements.append("VHT CAP Info, MAX MPDU Len %d, RX LDPC %d, TX STBC %d, RX STBC %d, Sounding Dim %d" % (
                            (tmpVhtCapInfo >> 0) & 3,
                            (tmpVhtCapInfo >> 4) & 1,
                            (tmpVhtCapInfo >> 7) & 1,
                            (tmpVhtCapInfo >> 8) & 7,
                            (tmpVhtCapInfo >> 16) & 7
                        ))
                    elif(tmpElement == MGMT_ELE.VHTOPS):
                        tmpMgmtElements.append("VHT Operation: To be added")
                    elif(tmpElement == MGMT_ELE.TXPOWER):
                        tmpMgmtElements.append("TX Power: Local Max Tx Pwr Count %d, Local Max Tx Pwr 20M %d" % (inbytes[elementIter] & 3, inbytes[elementIter+1]))
                    elif(tmpElement == MGMT_ELE.VENDOR):
                        tmpMgmtElements.append("Vendor: To be added")
                    else:
                        pass
                elementIter += tmpLen
            else:
                tmpElementByte = inbytes[elementIter]
                elementIter += 1
                tmpLen = inbytes[elementIter]
                elementIter += 1
                elementIter += tmpLen
        return tmpMgmtElements
    print("cloud mac80211header, mgmtParser input type error")
    return []


def pktParser(pkt):
    if(isinstance(pkt, (bytes, bytearray))):
        pktLen = len(pkt)
        procdLen = 0
        # fc
        if((procdLen + 2) <= pktLen):
            hdr_fc = frameControl(struct.unpack('<H', pkt[0:2])[0])
            hdr_fc.printInfo()
            procdLen += 2
        else:
            return
        # duration
        if((procdLen + 2) <= pktLen):
            hdr_duration = struct.unpack('<H', pkt[2:4])[0]
            procdLen += 2
            print("Packet duration %d us" % hdr_duration)
        else:
            return
        # check type
        if(hdr_fc.type == FC_TPYE.MGMT):
            if((procdLen + 18) <= pktLen):
                hdr_addRx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen:procdLen+6])
                hdr_addTx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen+6:procdLen+12])
                hdr_addDest = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen+12:procdLen+18])
                procdLen += 18
                print("Management to %s from %s dest %s" % (hdr_addRx, hdr_addTx, hdr_addDest))
            if((procdLen + 2) <= pktLen):
                hdr_seqCtrl= struct.unpack('<H', pkt[procdLen:procdLen+2])[0]
                procdLen += 2
                print("Management sequence control %d" % hdr_seqCtrl)
            if(hdr_fc.subType == FC_SUBTPYE_MGMT.BEACON):
                # Timestamp, Beacon Interval, Cap
                if((procdLen + 12) <= pktLen):
                    beacon_timestamp = struct.unpack('<Q', pkt[procdLen:procdLen+8])[0]
                    beacon_interval = struct.unpack('<H', pkt[procdLen+8:procdLen+10])[0]
                    beacon_cap = struct.unpack('<H', pkt[procdLen+10:procdLen+12])[0]
                    procdLen += 12
                    print("Management Beacon Timestamp %d, Interval %d, Cap %d" % (beacon_timestamp, beacon_interval, beacon_cap))
                # Elements
                if(procdLen <= pktLen):
                    beaconElements = mgmtElementParser(pkt[procdLen:])
                    for each in beaconElements:
                        print(each)
            elif(hdr_fc.subType == FC_SUBTPYE_MGMT.PROBEREQ):
                print("Management Probe Request")
            elif(hdr_fc.subType == FC_SUBTPYE_MGMT.PROBERES):
                print("Management Probe Response")
            else:
                print("Management subtype, not supported yet")
        elif(hdr_fc.type == FC_TPYE.CTRL):
            if(hdr_fc.subType == FC_SUBTPYE_CTRL.ACK):
                if((procdLen + 6) <= pktLen):
                    hdr_addRx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen:procdLen+6])
                    procdLen += 6
                    print("ACK to %s" % hdr_addRx)
            elif(hdr_fc.subType == FC_SUBTPYE_CTRL.BLOCKACK):
                if((procdLen + 12) <= pktLen):
                    hdr_addRx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen:procdLen+6])
                    hdr_addTx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen+6:procdLen+12])
                    procdLen += 12
                    print("BLOCK ACK to %s from %s" % (hdr_addRx, hdr_addTx))
                    # details to be added
            elif(hdr_fc.subType == FC_SUBTPYE_CTRL.RTS):
                if((procdLen + 12) <= pktLen):
                    hdr_addRx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen:procdLen+6])
                    hdr_addTx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen+6:procdLen+12])
                    procdLen += 12
                    print("RTS to %s from %s" % (hdr_addRx, hdr_addTx))
            elif(hdr_fc.subType == FC_SUBTPYE_CTRL.CTS):
                if((procdLen + 6) <= pktLen):
                    hdr_addRx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen:procdLen+6])
                    procdLen += 6
                    print("CTS to %s" % hdr_addRx)
            else:
                print("Control subtype, not supported yet")
        elif(hdr_fc.type == FC_TPYE.DATA):
            if((procdLen + 18) <= pktLen):
                hdr_addRx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen:procdLen+6])
                hdr_addTx = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen+6:procdLen+12])
                hdr_addDest = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[procdLen+12:procdLen+18])
                procdLen += 18
                print("Data to %s from %s dest %s" % (hdr_addRx, hdr_addTx, hdr_addDest))
            if((procdLen + 2) <= pktLen):
                hdr_seqCtrl= struct.unpack('<H', pkt[procdLen:procdLen+2])[0]
                procdLen += 2
                print("Data sequence control %d" % hdr_seqCtrl)
            if(hdr_fc.subType == FC_SUBTPYE_DATA.DATA):
                print("Data")
                if(procdLen <= pktLen):
                    print(pkt[procdLen:].hex())
            elif(hdr_fc.subType == FC_SUBTPYE_DATA.QOSDATA):
                print("QoS Data")
                if((procdLen+2) <= pktLen):
                    print("QoS control " + pkt[procdLen:procdLen+2].hex())
                    procdLen += 2
                if(procdLen <= pktLen):
                    print(pkt[procdLen:].hex())
            elif(hdr_fc.subType == FC_SUBTPYE_DATA.QOSNULL):
                print("QoS Null Data")
                if((procdLen+2) <= pktLen):
                    print("QoS control " + pkt[procdLen:procdLen+2].hex())
                    procdLen += 2
            else:
                print("Data subtype, not supported yet")
        else:
            print("cloud mac80211header, type not supported yet")



            
