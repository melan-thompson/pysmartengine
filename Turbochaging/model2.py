import sys

sys.path.append("../")
sys.path.append("../Properties of working fluids/")
from Cylinder import *

if __name__ == "__main__":
    import json

    with open("model2input.json", mode='r', encoding='UTF-8') as f:
        setting = json.load(f)

    from pandas import read_excel

    data = read_excel(setting["文件"], index_col="机型").loc[setting["机型"]]
    cyl = CylinderGeometry(data["bore,mm"] * 1.e-3, data["stroke,mm"] * 1.e-3, data["connecting rod length,mm"] * 1.e-3,
                           data["compression ratio"], data["number of cylinders"])

    Miller = MillerCycle(cyl, setting["环境压力,bar"]*1.e5, setting["环境温度,K"])

    n=setting["发动机转速,rpm"]
    # 压缩过程
    Miller.initCompressor(setting["压气机等熵效率"],setting["压气机压比"])

    print("压气机后压力{}bar".format(Miller.pk/1.e5))
    print("压气机后温度{}K".format(Miller.Tk))

    #进气过程
    Miller.intake(setting["进气门关闭角,BBDC"],setting["充量系数"])
    print("进气质量{}kg".format(Miller.data.table[4].data[-1]))
    print("进气流量{}kg/s".format(Miller.data.table[4].data[-1]/(30*4/n)))

    #压缩过程
    Miller.adiabaticompress()

    # 预混燃烧过程
    Miller.preBurn(setting["预混燃烧比例"],setting["过量空气系数"])
    # print("预混燃烧后压力{}bar".format(Miller.p4/1.e5))
    # print("预混燃烧后温度{}K".format(Miller.T4))

    #扩散燃烧比例
    Miller.disffBurn()
    print("燃烧结束角{}CA".format(Miller.data.table[0].data[-1]))

    # 膨胀过程
    Miller.expansion()

    # 排气过程
    Miller.supercriticalExhaust(setting["涡轮等熵效率"])

    Miller.subcritical()

    print("涡前压力{}bar".format(Miller.pit*((setting["环境压力,bar"]+setting["排气背压,bar"]))))
    print("涡前温度{}K".format(Miller.T7))

    Miller.analyze(plot=True,speed=setting["发动机转速,rpm"],index=setting["画图"])