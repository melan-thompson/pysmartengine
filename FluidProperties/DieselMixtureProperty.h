#ifndef DieselMixtureProperty
#define DieselMixtureProperty
#include "ArrayTable.h"

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
#endif // !AirProperty

