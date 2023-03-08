# Cloud Multi-User MIMO (CMU)

Introduction
------------
- The setup includes one AP (USRP B210 2x2) and two stations (USRP B200 1x1)
- The offical procedures in the IEEE 802.11 standard are:
    - 1. AP announces the channel sounding to stations.
    - 2. AP sends NDP.
    - 3. Stations estimate the channel.
    - 4. Stations compress the channel into feedback matrix V.
    - 5. Station 0 sends V to AP.
    - 6. AP polls station 1.
    - 7. Station 1 sends V to AP.
    - 8. AP generates beamforming steering matrix Q.
    - 9. AP apply Q to the spatial streams to send multi-user packets.


Version 1
-------------
- At the beginning, we try to use the full estimated channel on the station side, and we also assume that the CPU is not powerful enough to support 2 TX and 2 RX processing at the same time due to many USRP overflows **O** on the receiver side. In that case, we only let the AP to do the 2x2 signal file TX and the stations only need to do SISO RX. This is to test the gr-ieee80211 receiver, channel estimation, channel feedback and zero-forcing spatial mapping.
- The procedures change to
    - 1. AP sends NDP.
    - 2. Stations estimate the channel.
    - 3. Stations save channel into files
    - 4. AP gets channel files.
    - 5. AP generates beamforming steering matrix Q.
    - 6. AP apply Q to the spatial streams to send multi-user packets.
- This demo is in the folder **tools/cmu_v1** and **examples/cmu_v1**. Here are the steps:
    - 1. On the AP side, use **tools/pktGenExample.py** to generate 2x2 NDP signal file.
    - 2. On the station sides, run **tools/cmu_v1/cmu_sta0.py** and **tools/cmu_v1/cmu_sta1.py**.
    - 3. Find a clean channel.
    - 4. On the station sides, run **examples/cmu_v1/rxSta0.grc** and **examples/cmu_v1/rxSta1.grc** in GNU Radio, so far the stations are running, let's go to the AP side.
    - 5. AP side, use **examples/cmu_v1/txNdp.grc** in GNU Radio to send the NDP to stations.
    - 3. AP side, use **tools/cmu_v1/cmu_ap.py** to get the channel files and generate MU-MIMO packet.
    - 4. AP side, use **examples/cmu_v1/txMu.grc** in GNU Radio to send MU-MIMO packet.

### Notifications
- Choose a proper gain for USRP source and sink depending on the distance.
- For the **cmu_ap.py**, **cmu_sta0.py** and **cmu_sta1.py**, please change the information according to your own setting (like the IP, user name, password and path).
- The transmissions of the NDP and MU-MIMO packets could fail due to interference or noise. Try to resend if so.

Version 2
-------------
- The procedure is similar to the Version 1, but using the gr-ieee80211 transmitter to send the UDP and MU-MIMO packets. This is to test the gr-ieee80211 transmitter.
- This demo is in the folder **tools/cmu_v2** and **examples/cmu_v2** for the AP, while for the stations, still uses the files in **tools/cmu_v1** and **examples/cmu_v1**. Here are the steps:
    - 1. On the station sides, run **tools/cmu_v1/cmu_sta0.py** and **tools/cmu_v1/cmu_sta1.py**, the same as Version 1.
    - 2. Find a clean channel.
    - 3. On the station sides, run **examples/cmu_v1/rxSta0.grc** and **examples/cmu_v1/rxSta1.grc** in GNU Radio, the same as Version 1.
    - 4. On the AP side, run **examples/cmu_v2/txAp.grc**.
    - 5. On the AP side, run **tools/cmu_v2/cmu_ap.py**.
- The AP MAC first sends NDP packets, gets the channel files from the stations, updates the spatial mapping matrix to the GNU Radio transmitter and sends the MU-MIMO packet. This time, the AP only has MAC layer operations without PHY packet generations.

Version 3
-------------
- In this version, we no longer use the file to exchange channel information but set the wireless transceiver at both AP and station sides.
- This demo is in the folder **tools/cmu_v3** and **examples/cmu_v3** for the AP and stations. Here are the steps:
    - 1. On the station sides, run **tools/cmu_v3/cmu_sta0.py** and **tools/cmu_v3/cmu_sta1.py**.
    - 2. Find a clean channel.
    - 3. On the station sides, run **examples/cmu_v3/trxSta0.grc** and **examples/cmu_v3/trxSta1.grc** in GNU Radio.
    - 4. On the AP side, run **examples/cmu_v3/trxAp.grc**.
    - 5. On the AP side, run **tools/cmu_v3/cmu_ap.py**.
- For the stations side, everytime they receive the NDP to do the channel sounding, they put the channel info into a data packet and then send it to the AP. They just keep running for each NDP packet.
- For the AP side, it first sends an NDP packet and then wait for the data packet containing the channel info.
- The NDP will be retransmitted if the reception is timeout until it gets the channel info.
- The AP generates the Q matrix and sends MU-MIMO packets for multiple times.




# Beacon Transmission

Introduction
------------
- Generate a beacon packet which could be received by phones.

Steps
------------
- 1. Use the pktGenExample.py, the **SISO Legacy beacon example** part, generate a beacon packet signal bin file.
- 2. Use the txBeaconBin.grc in the **example/beacon** folder to transmit the bin file.
- 3. Pay attention on the file path of signal generation and file sink block in the grc.