import os
import struct
from matplotlib import pyplot as plt
import numpy as np
pyToolPath = os.path.dirname(__file__)
inSig = b''
addr = os.path.join(pyToolPath, "../tmp/sisorec.bin")



file_path = addr
chunk_size = 1024 * 1024  # 1MB

def myPlot(inSigComp):
    plt.plot(np.real(inSigComp))
    plt.plot(np.imag(inSigComp))
    plt.show()

def process_chunk(chunk):
    global inSig
    # Add your processing logic here
    print(f'Processing chunk of size: {len(chunk)} bytes')
    inSig+=chunk

count = 0
try:
    with open(file_path, 'rb') as file:
        while True:
            chunk = file.read(chunk_size)
            if not chunk:
                break
            process_chunk(chunk)
            count +=1
            if count == 10:
                break
            
except FileNotFoundError:
    print(f"File not found: {file_path}")
except IOError as e:
    print(f"IO error occurred: {e}")

print("len inSig",len(inSig))

sigSampNum = int(len(inSig) / 8)
sigComp = []
for i in range(0, sigSampNum):
    sigComp.append(complex(struct.unpack('f', inSig[i*8:i*8+4])[0], struct.unpack('f', inSig[i*8+4:i*8+8])[0]))


myPlot(sigComp)
# binF = open("recordchop", "wb")

# for i in range(0, sigSampNum):
#     binF.write(struct.pack("f", np.real(sigComp[i])))
#     binF.write(struct.pack("f", np.imag(sigComp[i])))
# binF.close()