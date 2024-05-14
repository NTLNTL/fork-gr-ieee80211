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

import mac80211header as m8h
import numpy as np
import struct
import socket
import binascii
import zlib

"""
    the simulation generates the IEEE 802.11a/g/n/ac MAC data packet and PHY data signal from transportation layer (UDP)
    ISO layer       |   protocol
    Transport       |   UDP
    Network         |   IPv4
    MAC             |   IEEE80211
    PHY             |   IEEE80211n
"""

def procCheckCrc32(inBytes):
    if(len(inBytes) > 4):
        inBytesPayload = inBytes[:-4]
        inBytesCrc32 = inBytes[-4:]
        tmpPayloadCrc32 = zlib.crc32(inBytesPayload)
        tmpTailCrc32, =struct.unpack('<L',inBytesCrc32)
        if(tmpPayloadCrc32 == tmpTailCrc32):
            return True
        else:
            return False
    else:
        print("cloud mac80211: crc32 input len error")

def procGenBitCrc8(bitsIn):
    c = [1] * 8
    for b in bitsIn:
        next_c = [0] * 8
        next_c[0] = b ^ c[7]
        next_c[1] = b ^ c[7] ^ c[0]
        next_c[2] = b ^ c[7] ^ c[1]
        next_c[3] = c[2]
        next_c[4] = c[3]
        next_c[5] = c[4]
        next_c[6] = c[5]
        next_c[7] = c[6]
        c = next_c
    return [1 - b for b in c[::-1]]

# udp generator, input: ip, port and payload
class udp():
    def __init__(self, sIp, dIp, sPort, dPort):
        self.fakeSourIp = sIp   # ip of network layer, only to generate checksum
        self.fakeDestIp = dIp
        self.sourPort = sPort   # ports
        self.destPort = dPort
        self.payloadBytes = b''
        self.protocol = socket.IPPROTO_UDP
        self.len = 8 
        self.checkSum = 0

    def __genCheckSum(self):
        self.checkSum = 0
        fakeSourIpBytes = socket.inet_aton(self.fakeSourIp)
        fakeDestIpBytes = socket.inet_aton(self.fakeDestIp)
        self.checkSum += fakeSourIpBytes[0] * 256 + fakeSourIpBytes[1]
        self.checkSum += fakeSourIpBytes[2] * 256 + fakeSourIpBytes[3]
        self.checkSum += fakeDestIpBytes[0] * 256 + fakeDestIpBytes[1]
        self.checkSum += fakeDestIpBytes[2] * 256 + fakeDestIpBytes[3]
        self.checkSum += self.protocol
        self.checkSum += self.len        # fake len = len
        self.checkSum += self.sourPort
        self.checkSum += self.destPort
        self.checkSum += self.len        # len
        for i in range(0, int(np.floor(len(self.payloadBytes) / 2))):
            self.checkSum += ((self.payloadBytes[i * 2]) * 256 + (self.payloadBytes[i * 2 + 1]))
        if (len(self.payloadBytes) > int(np.floor(len(self.payloadBytes) / 2)) * 2):
            self.checkSum += (self.payloadBytes[len(self.payloadBytes) - 1] * 256)
        while (self.checkSum > 65535):
            self.checkSum = self.checkSum % 65536 + int(np.floor(self.checkSum / 65536))
        self.checkSum = 65535 - self.checkSum
    def genPacket(self, payloadBytes):
        self.payloadBytes = payloadBytes
        self.len = len(payloadBytes) + 8
        self.__genCheckSum() #4332
        print("sour pack",struct.pack('>H',self.sourPort))
        return struct.pack('>HHHH',self.sourPort,self.destPort,self.len,self.checkSum)+self.payloadBytes

# ipv4 generator, takes UDP as payload, give id,
class ipv4():
    def __init__(self, id, ttl, sIp, dIp):
        self.ver = 4    # fixed, Version
        self.IHL = 5    # fixed, Internet Header Length, number of 32 bits
        self.DSCP = 0   # fixed
        self.ECN = 0    # fixed
        self.ID = id    # similar to sequence number
        self.flagReserved = 0   # fixed
        self.flagDF = 1      # fixed, do not frag
        self.flagMF = 0      # fixed, no more frag
        self.flag = (self.flagReserved << 2) + (self.flagDF << 1) + self.flagMF     # fixed
        self.fragOffset = 0     # no frag, so no frag offset
        self.TTL = ttl  # Time to live
        self.protocol = socket.IPPROTO_UDP  # fixed, protocol
        self.sourIp = sIp   # source ip
        self.destIp = dIp   # destination ip
        self.payload = b""
        self.len = 0
        self.checkSum = 0

    def __genCheckSum(self):
        self.checkSum = 0
        self.checkSum += (self.ver * (16 ** 3) + self.IHL * (16 ** 2) + self.DSCP * 4 + self.ECN)
        self.checkSum += self.len
        self.checkSum += self.ID
        self.checkSum += (self.flag * (2 ** 13) + self.fragOffset)
        self.checkSum += (self.TTL * 256 + self.protocol)
        sourIpBytes = socket.inet_aton(self.sourIp)
        destIpBytes = socket.inet_aton(self.destIp)
        self.checkSum += sourIpBytes[0] * 256 + sourIpBytes[1]
        self.checkSum += sourIpBytes[2] * 256 + sourIpBytes[3]
        self.checkSum += destIpBytes[0] * 256 + destIpBytes[1]
        self.checkSum += destIpBytes[2] * 256 + destIpBytes[3]
        while (self.checkSum > 65535):
            self.checkSum = self.checkSum % 65536 + int(np.floor(self.checkSum / 65536))
        self.checkSum = 65535 - self.checkSum

    def genPacket(self, payloadBytes):
        self.payload = payloadBytes     # payload, UDP packet
        self.len = self.IHL * 4 + len(payloadBytes) #20 +38 = 58B
        self.__genCheckSum()
        return struct.pack('>HHHHHH', (self.ver * (16 ** 3) + self.IHL * (16 ** 2) + self.DSCP * 4 + self.ECN), self.len, self.ID, (self.flag * (2 ** 13) + self.fragOffset), (self.TTL * 256 + self.protocol), self.checkSum) + socket.inet_aton(self.sourIp) + socket.inet_aton(self.destIp) + self.payload

# Logical Link Control, provide unified interface of data link layer
class llc():
    def __init__(self):
        self.SNAP_DSAP = 0xaa   # fixed, Source SAP
        self.SNAP_SSAP = 0xaa   # fixed, Destination SAP
        self.control = 0x03     # fixed
        self.RFC1024 = 0x000000 # fixed
        self.type = 0x0800  # 0x0800 for IP packet, 0x0806 for ARP

    def genPacket(self, payloadBytes):
        print("here",len(struct.pack('>BBB', self.SNAP_DSAP, self.SNAP_SSAP, self.control) + struct.pack('>L', self.RFC1024)[:3] + struct.pack('>H', self.type) + payloadBytes))
        return struct.pack('>BBB', self.SNAP_DSAP, self.SNAP_SSAP, self.control) + struct.pack('>L', self.RFC1024)[:3] + struct.pack('>H', self.type) + payloadBytes

class mac80211():
    """
    802.11 a & n mac frame
    | FC 2 | Duration 2 | ADDR1 6 | ADDR2 6 | ADDR3 6 | seq 2 | payload | FCS 4 |
    | FC 2 | Duration 2 | ADDR1 6 | ADDR2 6 | ADDR3 6 | seq 2 | QoS 2 | HT Control 4 | payload | FCS 4 |
    QoS field only appears in QoS packet: block ack, QoS data and so on
    HT Control field only appears in control wrapper
    QoS is added after a/b/g/, only used when the packet is QoS packet
    QoS field: 0-3: Traffic ID, used to show the priority of the packet, the other parts usually are 0
    """
    def __init__(self, type, subType, toDs, fromDs, retry, protected, addr1, addr2, addr3, seq):
        self.fc_protocol = 0        # fixed, frame control - protocol, 2 bits
        self.fc_type = type         # frame control - type, 2 bits
        self.fc_subType = subType   # frame control - sub type, 4 bits
        if(self.fc_subType == 8):
            print("cloud mac80211: QoS Data")
        elif(self.fc_subType == 0):
            print("cloud mac80211: Data")
        self.QoS = 0
        self.fc_toDs = toDs         # frame control - station to ap, 1 bit
        self.fc_fromDs = fromDs     # frame control - ap to station, 1 bit
        self.fc_frag = 0            # fixed, frame control - if more frag? 0: no more frags, 1 bit
        self.fc_retry = retry       # frame control - if retry? 1: retry, 0: not, 1 bit
        self.fc_pwr = 0    # fixed, frame control - not entering power saving, 1 bit
        self.fc_more = 0   # fixed, frame control - no data buffered, 1 bit
        self.fc_protected = protected   # frame control - if encrypted, 1 bit
        self.fc_order = 0  # fixed, frame control - no frag, no order, 1 bit
        self.duration = 0    # 16 bits, to be computed
        self.addr1 = addr1
        self.addr2 = addr2
        self.addr3 = addr3
        self.sc_frag = 0    # fixed, sequence control - frag number, no frag so to be 0
        self.sc_seq = seq     # sequence control - seq number

        self.sc = self.sc_frag + self.sc_seq << 4
        self.payloadBytes = b""
        self.fc = self.fc_protocol + (self.fc_type << 2) + (self.fc_subType << 4) + (self.fc_toDs << 8) + (self.fc_fromDs << 9) + (self.fc_frag << 10) + (self.fc_retry << 11) + (self.fc_pwr << 12) + (self.fc_more << 13) + (self.fc_protected << 14) + (self.fc_order << 15)
        self.eofPaddingSf = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.eofPaddingSf = self.eofPaddingSf + self.__macBitCrc8(self.eofPaddingSf) + [0, 1, 1, 1, 0, 0, 1, 0] #0x4e

    def __macBitCrc8(self, bitsIn):
        c = [1] * 8
        for b in bitsIn:
            next_c = [0] * 8
            next_c[0] = b ^ c[7]
            next_c[1] = b ^ c[7] ^ c[0]
            next_c[2] = b ^ c[7] ^ c[1]
            next_c[3] = c[2]
            next_c[4] = c[3]
            next_c[5] = c[4]
            next_c[6] = c[5]
            next_c[7] = c[6]
            c = next_c
        return [1 - b for b in c[::-1]]

    def __genDuration(self):
        # manually set
        self.duration = 110  # used in sniffed packet

    def genPacket(self, payload):
        if(isinstance(payload, (bytes, bytearray))):
            self.payloadBytes = payload
            self.__genDuration()
            tmpPacket = struct.pack('<H', self.fc) + struct.pack('<H', self.duration)
            tmpPacket += binascii.unhexlify("".join(self.addr1.split(":")))
            tmpPacket += binascii.unhexlify("".join(self.addr2.split(":")))
            tmpPacket += binascii.unhexlify("".join(self.addr1.split(":")))
            tmpPacket += struct.pack('<H', self.sc)
            if(self.fc_subType == 8):
                tmpPacket += struct.pack('<H', self.QoS)   # only added when QoS packet
            tmpPacket += self.payloadBytes
            tmpPacket += struct.pack('<L',zlib.crc32(tmpPacket))

            print("cloud mac80211, gen pkt mac mpdu length: %d" % len(tmpPacket))
            return tmpPacket
        else:
            print("cloud mac80211, gen pkt input type error")
            return b""
    
    def genCtrlVhtNdpAnnouncement(self, rxAddr, txAddr, token, staAids, staFbType, staNc):
        if(isinstance(rxAddr, str) and isinstance(txAddr, str) and isinstance(token, int) and isinstance(staAids, list) and isinstance(staFbType, list) and isinstance(staNc, list)):
            if(token >= 0 and token <= 63):
                tmpFc = 0 + (1 << 2) + (5 << 4) + (0 << 15)     # fc flag part all zero, type control, subtype vht ndp announcement
                tmpDialogToken = (token << 2) + 0               # 2 bits reserved
                tmpPacket = struct.pack('<H', tmpFc) + struct.pack('<H', 340)
                tmpPacket += binascii.unhexlify("".join(rxAddr.split(":")))
                tmpPacket += binascii.unhexlify("".join(txAddr.split(":")))
                tmpPacket += struct.pack('<B', tmpDialogToken)
                nSta = min([len(staAids), len(staFbType), len(staNc)])
                for i in range(0, nSta):
                    tmpStaInfo = staAids[i] + (staFbType[i] << 12)
                    if(staFbType[i]):
                        tmpStaInfo + ((staNc[i] - 1) << 13)
                    else:
                        tmpStaInfo + (0 << 13)      # for SU this is reserved
                    tmpPacket += struct.pack('<H', tmpStaInfo)
                tmpPacket += struct.pack('<L',zlib.crc32(tmpPacket))
                return tmpPacket
        print("cloud mac80211, vht ndp announcement input params error")
        return b""
    
    def genCtrlBfReportPoll(self, rxAddr, txAddr, fbSegList):
        if(isinstance(rxAddr, str) and isinstance(txAddr, str) and isinstance(fbSegList, list)):
            if(len(fbSegList) < 9 and min(fbSegList) >= 0 and max(fbSegList) <= 7):
                tmpFc = 0 + (1 << 2) + (4 << 4) + (0 << 15)     # fc flag part all zero, type 1 control, subtype 4 bf report poll
                tmpPacket = struct.pack('<H', tmpFc) + struct.pack('<H', 110)
                tmpPacket += binascii.unhexlify("".join(rxAddr.split(":")))
                tmpPacket += binascii.unhexlify("".join(txAddr.split(":")))
                tmpFbSegBitmap = 0
                for each in fbSegList:
                    tmpFbSegBitmap |= (1 << each)
                tmpPacket += struct.pack('<B', tmpFbSegBitmap)
                tmpPacket += struct.pack('<L',zlib.crc32(tmpPacket))
                return tmpPacket
        print("cloud mac80211, beamforming report poll input params error")
        return b""

    # ieee80211-2020 section 9.4.1.11
    def genMgmtActNoAck(self, dsAddr, txAddr, bssid, seq, category, details):
        if(isinstance(dsAddr, str) and isinstance(txAddr, str) and isinstance(bssid, str) and isinstance(details, (bytes, bytearray))):
            tmpFc = 0 + (0 << 2) + (14 << 4) + (0 << 15)     # fc flag part all zero, type 1 control, subtype 5 vht ndp announcement
            tmpPacket = struct.pack('<H', tmpFc) + struct.pack('<H', 32)
            tmpPacket += binascii.unhexlify("".join(dsAddr.split(":")))
            tmpPacket += binascii.unhexlify("".join(txAddr.split(":")))
            tmpPacket += binascii.unhexlify("".join(bssid.split(":")))
            tmpPacket += struct.pack('<H', seq)
            tmpPacket += struct.pack('<B', category)
            tmpPacket += details
            tmpPacket += struct.pack('<L',zlib.crc32(tmpPacket))
            return tmpPacket
        print("cloud mac80211, management action no ack input params error")
        return b""
    
    def mgmtActNoAckParser(self, pkt):
        if(isinstance(pkt, (bytes, bytearray))):
            if(len(pkt) > 25):
                hdr_fc = m8h.frameControl(struct.unpack('<H', pkt[0:2])[0])
                hdr_dsAddr = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[4:10])
                hdr_txAddr = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[10:16])
                hdr_bssid = "%x:%x:%x:%x:%x:%x" % struct.unpack("BBBBBB",pkt[16:22])
                mgmtActCategory = m8h.MGMT_ACT_CAT(int(pkt[24]))
                return mgmtActCategory, pkt[25:]
        return 0, b""
            

def genAmpduHT(payloads):
    if(isinstance(payloads, list)):
        tmpAmpduPkt = b""
        for payloadIter in range(0, len(payloads)):
            delimiterMpduLen = len(payloads[payloadIter])
            if(delimiterMpduLen < 1 or delimiterMpduLen > 4095):
                print("cloud mac80211, gen ampdu ht: packet %d len %d error" % (payloadIter, delimiterMpduLen))
                return b""
            delimiterMpduLenBits = []
            for i in range(0, 12):      # HT only 12 bits for len
                delimiterMpduLenBits.append((delimiterMpduLen >> i) & (1))
            delimiterBits = [0, 0, 0, 0] + delimiterMpduLenBits
            delimiterBits = delimiterBits + procGenBitCrc8(delimiterBits)
            for i in range(0, 8):
                delimiterBits.append((0x4e >> i) & (1))
            tmpDelimiterBytes = b""
            for i in range(0, 4):
                tmpByte = 0
                for j in range(0, 8):
                    tmpByte = tmpByte + delimiterBits[i * 8 + j] * (2 ** j)
                tmpDelimiterBytes += bytearray([tmpByte])
            tmpPacket = tmpDelimiterBytes + payloads[payloadIter]
            if(payloadIter < (len(payloads)-1)):  # pad if not last one
                nBytePadding = int(np.ceil(len(tmpPacket) / 4) * 4 - len(tmpPacket))
                tmpPacket += b'\x00' * nBytePadding
            tmpAmpduPkt += tmpPacket
        return tmpAmpduPkt
    else:
        print("cloud mac80211, gen ampdu ht: input type error")
        return b""

def genAmpduVHT(payloads):
    if(isinstance(payloads, list)):
        tmpAmpduPkt = b""
        for eachPayload in payloads:
            delimiterMpduLen = len(eachPayload)
            delimiterMpduLenBits = []
            delimiterEof = 0
            if(len(payloads) == 1):
                delimiterEof = 1
            delimiterReserved = 0
            for i in range(0, 14):
                delimiterMpduLenBits.append((delimiterMpduLen >> i) & (1))
            delimiterBits = [delimiterEof] + [delimiterReserved] + delimiterMpduLenBits[12:14] + delimiterMpduLenBits[0:12]
            delimiterBits = delimiterBits + procGenBitCrc8(delimiterBits)
            for i in range(0, 8):
                delimiterBits.append((0x4e >> i) & (1))
            tmpDelimiterBytes = b""
            for i in range(0, 4):
                tmpByte = 0
                for j in range(0, 8):
                    tmpByte = tmpByte + delimiterBits[i * 8 + j] * (2 ** j)
                tmpDelimiterBytes += bytearray([tmpByte])
            tmpPacket = tmpDelimiterBytes + eachPayload
            # each packet is padded
            nBytePadding = int(np.ceil(len(tmpPacket) / 4) * 4 - len(tmpPacket))
            tmpPacket += b'\x00' * nBytePadding
            tmpAmpduPkt += tmpPacket
        return tmpAmpduPkt
    else:
        print("cloud mac80211, gen ampdu vht: input type error")
        return b""



class RxMac80211():
    def __init__(self,inBitArray,debug):
        self.dbFlag = debug
        self.procIndex = 0 
        self.fc_subType = self.Bi2Dec(inBitArray[4:8])


        self.macLen = 28 # bytes 

        if self.fc_subType == 0:
            self.macLen = 24
        elif self.fc_subType == 8:     #sub type, 8 = QoS Data, 0 = Data
            self.macLen = 26
        else:
            print("ERROR:unsupport sub type [%d]"%(self.fc_subType))
            return
        self.macBits= inBitArray
        self.llcBits = self.macBits[self.macLen*8:]
        self.IpBits = self.llcBits[8*8:]
        self.udpBits = self.IpBits[20*8:]

        #output
        self.udpPayload = ""
    def rxMacStepsList(self):
        self.__procRxMacMac()
        self.__procRxLLCLLC()
        self.__procRxIPIP()
        self.__procRxUDP()
        return self.udpPayload
    def __procRxMacMac(self):
        fc_protocol = self.Bi2Dec(self.macBits[0:2])
        fc_type = self.Bi2Dec(self.macBits[2:4])
        fc_subType = self.fc_subType
        fc_toDs = self.Bi2Dec(self.macBits[8])
        fc_fromDs = self.Bi2Dec(self.macBits[9])
        fc_frag = self.Bi2Dec(self.macBits[10])
        fc_retry = self.Bi2Dec(self.macBits[11])
        fc_pwr = self.Bi2Dec(self.macBits[12])
        fc_more = self.Bi2Dec(self.macBits[13])
        fc_protected = self.Bi2Dec(self.macBits[14])
        fc_order = self.Bi2Dec(self.macBits[15])
        duration = self.Bi2Dec(self.macBits[16:32])
        addr1 = self.Bi2hexSting(self.macBits[4*8:10*8])
        addr2 = self.Bi2hexSting(self.macBits[10*8:16*8])
        addr3 = self.Bi2hexSting(self.macBits[16*8:22*8])
        sc_frag = self.Bi2Dec(self.macBits[22*8:22*8+4])
        sc_seq = int(self.Bi2Dec(self.macBits[22*8:24*8])>>4)

        if self.fc_subType == 0:
            if self.dbFlag:print("|----------------------------------------------------MAC HEADER----------------------------------------------------|")
            if self.dbFlag:print("Frame Control: |Protocol Ver - %d|Type - %d|Subtype - %d|" %(fc_protocol,fc_type,fc_subType))
            if self.dbFlag:print("MAC Frame: |Duration - %d|Addr 1 - %s|Addr 2 - %s|Addr 3 - %s|sc_frag - %d|sc_seq - %d|"%(duration,addr1,addr2,addr3,sc_frag,sc_seq))
        elif self.fc_subType == 8:
            ht_control = None
            if self.dbFlag:print()
    def __procRxLLCLLC(self):
        if self.dbFlag:print("|----------------------------------------------------LLC HEADER----------------------------------------------------|")
        # print(self.llcBits)
        tmpSSAP = self.llcBits[0:8]
        print("tmpSSAP",self.Bi2hexSting(tmpSSAP))
    def __procRxIPIP(self):
        if self.dbFlag:print("|----------------------------------------------------IP HEADER-----------------------------------------------------|")
        # print(self.IpBits)
    def __procRxUDP(self):

        udpHeader = struct.unpack(">HHHH", self.Bi2Byte(self.udpBits[0:8*8]))
        sourPort = udpHeader[0]
        destPort = udpHeader[1]
        udpLen = udpHeader[2]
        udpCheckSum = udpHeader[3]


        tmpPayloadbit = self.udpBits[8*8:udpLen*8]
        tmpB = ""
        for i in range(int(len(tmpPayloadbit)/8)):
            tmp = tmpPayloadbit[i*8:i*8+8]
            tmpB+=chr(self.Bi2Dec(tmp))
        if self.dbFlag:print("|----------------------------------------------------UDP HEADER----------------------------------------------------|")
        # print(self.udpBits)
        self.udpPayload = tmpB
        # print("udp Payload:",tmpB)

    def Bi2Dec(self,inBi):
        tmpSum = 0
        if isinstance(inBi, list):
            for i in range(len(inBi)-1,-1,-1):
                if inBi[i] == 1:
                    tmpSum+= (1<<i)
        else:
            tmpSum = inBi
        return tmpSum
    def Bi2hexSting(self,inBi):
        tmpOut = ""
        tmpIpv6 = ""
        for i in range(int(len(inBi)/8)): #six bytes for addr field 
            tmp = hex(self.Bi2Dec(inBi[i*8:i*8+8]))[2:].rjust(2,'0')
            tmpOut += tmp
        tmpIpv6 = ":".join([tmpOut[i:i+2] for i in range(0, len(tmpOut),2)])
        return tmpIpv6
    def Bi2Byte(self,inBi): #the input of unpack is byte, the ouptut of this function is dec number
        tmpOut = b""
        if isinstance(inBi, list):
            
            for i in range(int(len(inBi)/8)):
                tmp = self.Bi2Dec(inBi[i*8:i*8+8])
                tmpOut += tmp.to_bytes(1,byteorder='big')
            return tmpOut
            # return (tmp.to_bytes(tmp,byteorder='big'))
            
        else:
            return "ERROR: in Bi2Byte2Dec"
if __name__ == "__main__":
    mac80211Ins = mac80211(2,  # type
                                     8,  # sub type, 8 = QoS Data
                                     1,  # to DS, station to AP
                                     0,  # from DS
                                     0,  # retry
                                     0,  # protected
                                     'f4:69:d5:80:0f:a0',  # dest add
                                     '00:c0:ca:b1:5b:e1',  # sour add
                                     'f4:69:d5:80:0f:a0',  # recv add
                                     2704)  # sequence
    print(mac80211Ins.genCtrlVhtNdpAnnouncement("6e:1b:72:2a:1c:b8", "00:27:e3:9d:e8:9c", 23, [100], [0], [0]).hex())
    