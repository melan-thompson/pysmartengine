# 项目介绍
本项目为内燃机性能仿真常用代码，用python语言开发，使用python编写的目的是能够使用最新的算法和框架。因为pytorch能够用GPU加速矩阵的运算，后期将使用pytorch框架构建一维非定常流动以及缸内CFD计算模块。本项目搭建的数据库均使用SQLite，数据交互使用json。

# 本项目的模块包括

- 工质热物性计算
- 气缸模块
  - 示功图的分析
  - 缸内工作过程计算
- 管道模块
  - 一维非定常流动模型
  - 管道快速计算QPM模型
  - MOR降阶模型
- 气阀模块
- 涡轮增压器模块
- 发动机性能模块
  - 发动机的平均值模型(Mean Value Engine Models)
  - 发动机的容积法模型(Filling and Emptying models)
  - 一维计算模型(1D engine models)
  - 一维快速计算模型
  - 一三维耦合模型



# 附录 常用的内燃机性能仿真软件介绍

[意大利OpenWAM](http://openwam.webs.upv.es/docs/) OpenWAM是一个开源的发动机仿真软件，[openWAM源代码](https://github.com/CMT-UPV/OpenWAM.git)

[俄罗斯DIESEL-RK](https://diesel-rk.bmstu.ru/Eng/index.php) DIESEL-RK is a full cycle thermodynamic engine simulation software. One is designed for simulating and optimizing working processes of **two**- and **four-stroke** internal combustion engines with all types of boosting.

[Ansys Chemkin](https://www.ansys.com/products/fluids/ansys-chemkin-pro) Ansys Chemkin-Pro is a chemical kinetics simulator that models idealized reacting flows and provides insight into results before production testing.

[AVL FIRE](https://www.avl.com/fire) 发动机三维流场仿真模拟

[AVL Boost](https://www.avl.com/boost/)



## 一维性能计算技术的发展

1964年以前，以热力循环分析为主。

1964年，英国曼彻斯特理工学院(UMIST)Renolds Benson教授首次运用特征线法对发动机进排气系统内的压力波动现象求解。

1965年，英国贝尔法斯特(Belfast)大学Golden Blair教授首次将一维特征线法应用于二冲程柴油机一维换气过程的模拟优化计算。

1993年，英国里卡多（Ricardo)咨询公司在美国推出Wave模拟软件，首次引入图形式用户界面（GUI），大大提高了发动机一维性能模拟软件的使用方便性。

1996年伽玛技术公司（Gamma Technologies,Inc.) 推出发动机一维性能模拟软件GT-Power。

1998年左右AVL推出Boost。

## 一维实时仿真

