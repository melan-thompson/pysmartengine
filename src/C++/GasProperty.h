#ifndef GASPROPERTY
#define GASPROPERTY
#include <iostream>
#include <cmath>
#include "ArrayTable.h"

namespace justi
{
	double Rg(double AFAK);//J.keenan & J.Kaye formula, used for ideal gas
	

	double u_Justi(double T, double AFAK = 1.e8, double Tref = 273.15);


	double h_Justi(double T, double AFAK = 1.e8, double Tref = 273.15);
	

	double cv_Justi(double T, double AFAK = 1.e8, double Tref = 273.15);
	

	double M(double AFAK);//J.keenan & J.Kaye formula, used for ideal gas
	

	double k_Justi(double T, double AFAK = 1.e8);
	

	double cp_Justi(double T, double AFAK = 1.e8, double Tref = 273.15);
	

	double cv_mean_Justi(double T_begin, double T_end, double AFAK = 1.e8, double step = 0.1);

	double cp_mean_Justi(double T_begin, double T_end, double AFAK = 1.e8, double step = 0.1);

	double k_mean_Justi(double T_begin, double T_end, double AFAK = 1.e8, double step = 0.1);

}

namespace schmidt//calculate cv using F,Schmidt's fomula
{
	double Mcv_pure_air(double T);//T(k),M*cv(kJ/(kmol*K))

	double Mcv_pure_exhaust_gas(double T);//T(k),M*cv(kJ/(kmol*K))
}

namespace FZacharias
{
	double Z(double p, double T, double AFAK);//采用F.Zacharias公式计算压缩因子，已验证其正确性
		//p(Pa),T(k),AFAK(1),Z(1)
	

	double Rg(double p, double T, double AFAK);//used for high pressure, F.Zacharias formula
		//p(Pa),T(K),AFAK(1),Rg(J/(kg*K))
	
}


class DieselMixture
{
public:
	double ma, r,xr;//r=mr/mc; 需要初始化，初始化后不发生变化
	double L0;//常量
	double gf;
	ArrayTable property;

	DieselMixture() {};

	///////默认的初始化方法，给出空气的量，残余废气系数
	DieselMixture(double ma, double r, double gf);

	bool initWithMcAFA(double _mc, double _afa, double _gf);

	bool initWithMcMa(double _mc, double _ma, double _gf);

	
	double AFAK(double x);

	double M_total(double x);

	double M_exhaust(double x);
	

	double M_air(double x);
	

	double Xk(double x);
	

	double cv(double x, double T);
	

	double cp(double x, double T);
	

	double U(double x, double T);
	

	double Rg(double x);
	

	double k(double x, double T);
	
	
private:

	bool _initPrepertyTable();
};

#endif // !GASPROPERTY

