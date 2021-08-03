#ifndef GASPROPERTY
#define GASPROPERTY
#include <iostream>
#include <cmath>

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

#endif // !GASPROPERTY

