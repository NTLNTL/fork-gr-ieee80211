
import phy80211header as p8h



def procRxNonDataSC(inQam): #remove non-data subcarrier
    if(isinstance(inQam, list) and len(inQam) in [64, 128, 256]):
        if(len(inQam) in [64,128,256]):
            return inQam[6:-5]
        else: #need verify 
            print("cloud rx phy80211 header need verified")
    else:
        print("cloud rx phy80211 header, procRxNonDataSC: input length error %d" % (len(inQam)))
        return []




if __name__ == "__main__":

    pass
