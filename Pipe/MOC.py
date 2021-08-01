def RiemannConst(a,v,k=1.4):
    lambdaa=a+(k-1)/2*v
    beta=a-(k-1)/2*v
    return lambdaa,beta

def Characteristics(lambdaa:float,beta:float,k:float=1.4):
    """
    返回三条特征线的方向
    :param lambdaa:dx/dt=v+a特征线
    :param beta:dx/dt=v-a特征线
    :param k:等熵指数
    :return:
    """
    temp1=(k+1)/2/(k-1)
    temp2=(3-k)/2/(k-1)
    Clambda=temp1*lambdaa-temp2*beta
    Cbeta=temp2*lambdaa-temp1*beta
    Cm=(lambdaa-beta)/(k-1)


class Node:
    def __init__(self,u,p,T,k=1.4,Rg=287):
        self.k=k
        self.Rg=Rg
        self.u=u

        self.lambdaa=None #左特征线黎曼不变量
        self.beta=None #右特征线黎曼不变量
        self.rho=p/Rg/T

        from math import sqrt
        self.a=sqrt(k*p/self.rho)


    def dlambda(self):
        pass

    def _detaRa(self,A,DADX):
        """
        管截面变化对黎曼不变量的影响
        :param A: 管截面
        :param DADX: 管截面的变化率
        :return:
        """
        return -(self.k-1)/2.*self.a*self.u/A/DADX

    def _Rh(self,q):
        """
        管路摩擦对黎曼不变量的影响
        :param q:单位质量的传热
        :return:
        """
        return pow(self.k-1,2)/2/self.a*q