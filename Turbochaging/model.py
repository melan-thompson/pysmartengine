
import sys
sys.path.append("../")
sys.path.append("../Properties of working fluids/")
from Cylinder import *
from Compressor import *
from GasProperty import Rg,k_Justi

if __name__=="__main__":
    import json

    with open("input.json", mode='r', encoding='UTF-8') as f:
        setting = json.load(f)

    from pandas import read_excel

    data = read_excel(setting["文件"], index_col="机型").loc[setting["机型"]]
    cyl = CylinderGeometry(data["bore,mm"] * 1.e-3, data["stroke,mm"] * 1.e-3, data["connecting rod length,mm"] * 1.e-3,
                           data["compression ratio"], data["number of cylinders"])

    # cyl.plotVolume()
    n=setting["发动机转速,rpm"]

    # 计算平均有效压力
    pme=BMEP(n,setting["发动机功率,kW"]*1.e3,cyl.totalDisplacedVolume())

    #计算活塞平均速度
    cm=cyl.stroke*n/30
    print("活塞平均速度{}m/s".format(cm))

    #计算机械效率
    etam=MeEff(cyl.bore,cm,pme)
    print("机械效率{}".format(etam))

    #计算压气机后压力和温度
    Pk=setting["环境压力,bar"]*setting["压气机压比"]*1.e5
    Tk=TAfterCompressor(setting["环境温度,K"],setting["压气机压比"],setting["压气机等熵效率"])
    Ts=Tk-setting["中冷器温降,K"]

    print("压气机后温度{}".format(Tk))
    print("进气管温度{}".format(Ts))


    # 计算涡轮前热利用系数
    R=thermalUsageIndex(setting["燃油消耗率,g/(kW*h)"]*etam,setting["过量空气系数"],cm,cyl.bore,cyl.stroke,n,Pk,Ts)
    print("涡轮前热利用系数{}".format(R))

    # 计算涡前温度
    Tt=exhaustTemperature(setting["燃油消耗率,g/(kW*h)"],etam,setting["过量空气系数"],setting["扫气系数"],Ts,R)
    print("涡前温度{}K".format(Tt))

    #涡轮后压力
    Pt0=(setting["环境压力,bar"]+setting["排气背压,bar"])*1.e5

    #计算进气量
    ma=massAir(cyl.totalDisplacedVolume(),Pk,Ts,setting["充量系数"],setting["扫气系数"])/(30*4/n)
    print("进气流量:{}".format(ma))

    #计算涡轮膨胀比
    from Valve import flowUnitArea

    def fun(pt):
        Rge=Rg(setting["过量空气系数"])
        k=k_Justi(Tt,setting["过量空气系数"])
        temp=setting["涡轮当量流通面积,m^2"]*flowUnitArea(pt,Tt,Rge,k,Pt0)
        print("涡前压力{},流量{}".format(pt,temp))
        return (temp-ma)/ma

    pt=2*Pt0
    h=1.
    while abs(fun(pt))>1.e-5:
        pt-=fun(pt)/((fun(pt+h)-fun(pt-h))/2./h)
        # print((fun(pt)))

    print("迭代完成，涡前压力{}bar".format(pt/1.e5))

    # 计算压气机压比
    pk=piK(Tt,pt/Pt0,setting["涡轮等熵效率"],setting["压气机等熵效率"],setting["环境温度,K"],setting["过量空气系数"]*14.3,ma)
    print(pk)



