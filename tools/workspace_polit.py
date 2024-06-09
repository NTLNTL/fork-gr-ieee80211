import numpy as np
import cmath 
import math
from matplotlib import pyplot as plt

def myConstellationPlot(inSig):
        x = []
        y = [] 
        absArray = []
        for i in range(len(inSig)):
            x.append(np.real(inSig[i]))
            y.append(np.imag(inSig[i]))
            absArray.append(abs(inSig[i]))
        # myPlot(absArray)
        plt.title("ConstellationPlot in self")
        plt.xlim([-1.1, 1.1])
        plt.ylim([-1.1, 1.1])
        plt.scatter(x,y)
        plt.show()

std_polit_p = [1+0j,1+0j,1+0j,-1+0j]

rcv_std_polit = [1.01+0.01j,1.11+0.11j,1.21+0.21j,-1.1+0.01j]


rcv_sample = 2.0+0.01j
diff1 = 0
diff2 = 0
for i in range(4):
    diff1 += std_polit_p[i] * np.conj(rcv_std_polit[i])
    diff2 += std_polit_p[i] / rcv_std_polit[i]
print(diff1,diff2)


print("new ------")
inQam = [(0.9420406688567533-0.4863904751472493j), (-0.8815589184953447+0.4215802580712356j), (0.9305685727498846-0.4205618906972427j), (0.9311463842420171-0.4193080612413801j), (-0.8807572884226389+0.42839270069715185j), (0.9688115614742474-0.46157188019080253j), (-0.8419513490035313+0.3911870731166867j), (0.9755044079860196-0.4483779493583563j), (-0.8385244017926303+0.40464457978588736j), (-0.8872018248521603+0.4401375619355113j), (-0.8887098166054521+0.44185101488787276j), (0.9807992501777204-0.41935701184540025j), (-0.8938385655177148+0.4454665947625716j), (-0.8936877800729717+0.4469909110243424j), (-0.8989318331506521+0.44683825624635104j), (-0.8417079770263467+0.45115018147453j), (0.9091540335689517-0.3984209248518937j), (0.9655530287828972-0.3781678478824294j), (0.9604324117509039-0.3724236600215625j), (0.9023819153615393-0.39817209419937954j), (0.948656522185681-0.36324143552142774j), (-0.8667612719672377+0.48231746918519264j), (0.9364425554488913-0.35615250237983953j), (-0.9200342110153108+0.4443845731193839j), (0.9228414040121891-0.35229346115543125j), (0.9151001630881851-0.3512053571467243j), (-0.9279111152344552+0.4364502312139099j), (0.8843006550936575-0.4098527338198162j), (0.8833446377145183-0.41121691064346655j), (0.8821246966919017-0.35615643830444565j), (0.8746823636032331-0.35826155922204084j), (0.8813713031438597-0.41854933458129445j), (0.8627823214048455-0.3663899419006512j), (0.8802383126286316-0.42424235665791044j), (0.8545804427673281-0.375215164178125j), (-0.9640958972484153+0.46768623042025054j), (-0.9685379995410073+0.4623845244782403j), (-0.9741364658313+0.4559508938460317j), (-0.9755668241175587+0.44919495779425545j), (-0.9778588431868536+0.4418809768626897j), (0.8374247549507995-0.4121342118925489j), (0.8350996183335291-0.41786839718481583j), (0.8915401712028247-0.4442761537448913j), (-0.9802459082219838+0.4121006840824481j), (-0.9170258395179687+0.400527827978371j), (0.8986923410172847-0.4474960132069459j), (-0.9111170468902059+0.3987760484188368j), (0.9034274683998431-0.448860778817892j), (0.8478681208675702-0.4629731901240341j), (0.8515566762653748-0.4689945698566545j), (0.8566896462230309-0.4735967868759957j), (0.8618672107413746-0.47781792173010756j)]
l = [(0.9688115614742474-0.46157188019080253j), (0.9023819153615393-0.39817209419937954j) ,(0.8627823214048455-0.3663899419006512j) ,(-0.9111170468902059+0.3987760484188368j)]
sum = 0 
std = [1,1,1,-1]
for i in range(4):
     sum+= (l[i]*np.conj(std[i]))
abssum = np.abs(sum)
diff = sum/abssum
print(sum,diff)
out = [each / diff  for each in l]
print(out)


#  diff = std_polit/rcv_std_polit
# recover  = diff * rcv_sample



# two_diff = b * np.conj(a)
# polar_two = cmath.polar(two_diff)
# print(two_diff,polar_two)
# # d = b/polar_two 
# print("diff and two_diff", diff ,two_diff)


# theta = math.atan2(two_diff.imag, two_diff.real) * 180/math.pi
# print(math.atan2(two_diff.imag, two_diff.real),theta)
# myConstellationPlot(l)