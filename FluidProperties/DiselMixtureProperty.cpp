#include "GasProperty/DieselMixtureProperty.h"
#include "GasProperty/GasProperty.h"


DieselMixture::DieselMixture(double ma, double r, double gf) :ma(ma), r(r), gf(gf)
{
	L0 = 14.3;
	if (ma < L0 * gf)
	{
		std::cerr << "Charged air is too few to burn all the fuel!!\n";
		std::cerr << "There are must at least " << L0 * gf << "(kg) air!!\n";
		system("pause");
	}

	if (r < 0 || r>1)
	{
		std::cerr << "residual gas fraction can not less than 0 and greater than 1!!\n";
		system("pause");
	}

	xr = ma * r / (1 - r) / (1 + L0) / gf;
	_initPrepertyTable();

	//AFAK = (L0 *gf* xr + ma) / (L0 *gf* xr);
};

bool DieselMixture::initWithMcAFA(double _mc, double _afa, double _gf)
{
	L0 = 14.3;
	gf = _gf;
	if (_afa < 1)
	{
		std::cerr << "Excess air fuel ratio must be greater than or equal to 1\n";
		system("pause");
	}
	double mr = (1 + L0) / (L0 * _afa + 1) * _mc;
	ma = L0 * (_afa - 1) / (L0 * _afa + 1) * _mc;
	if (ma < L0 * gf)
	{
		std::cerr << "Charged air is too few to burn all the fuel!!\n";
		std::cerr << "There are must at least " << L0 * gf << "(kg) air!!\n";
		system("pause");
	}
	r = mr / (_mc);
	xr = ma * r / (1 - r) / (1 + L0) / gf;

	_initPrepertyTable();
	return true;
}

bool DieselMixture::initWithMcMa(double _mc, double _ma, double _gf)
{
	L0 = 14.3;
	ma = _ma; gf = _gf;
	if (ma < L0 * gf)
	{
		std::cerr << "Charged air is too few to burn all the fuel!!\n";
		std::cerr << "There are must at least " << L0 * gf << "(kg) air!!\n";
		system("pause");
	}
	double mr = _mc - _ma;
	r = mr / _mc;
	xr = ma * r / (1 - r) / (1 + L0) / gf;
	_initPrepertyTable();
	return true;
}


double DieselMixture::AFAK(double x)
{
	if (x < 0)
	{
		x = 0;
	}
	/*else if (x > 1)
	{
		x = 1;
	}*/
	if ((xr + x) < 1.e-8)
	{
		return 1.e8;
	}
	return (L0 * gf * xr + ma) / (L0 * gf * (xr + x));
}

double DieselMixture::M_total(double x)
{
	if (x < 0)
	{
		x = 0;
	}
	/*else if (x > 1)
	{
		x = 1;
	}*/
	if ((xr + x) < 1.e-8)
	{
		return ma;
	}
	return gf * (xr + x) * (1 + AFAK(x) * L0);
}

double DieselMixture::M_exhaust(double x)
{
	return M_total(x) * (1 + L0) / (L0 * AFAK(x) + 1);
}

double DieselMixture::M_air(double x)
{
	return M_total(x) - M_exhaust(x);
}

double DieselMixture::Xk(double x)
{
	return M_total(x) / (1 + L0 * AFAK(x)) / gf;
	//Xk = M / (1 + Inj.L0 * AFA) / Inj.MF;
}

double DieselMixture::cv(double x, double T)
{
	using namespace justi;
	return cv_Justi(T, AFAK(x));
}

double DieselMixture::cp(double x, double T)
{
	using namespace justi;
	return cp_Justi(T, AFAK(x));//k(T, AFAK(x)) * cv(x, T);
}

double DieselMixture::U(double x, double T)
{
	using namespace justi;
	return M_total(x) * u_Justi(T, AFAK(x));
}

double DieselMixture::Rg(double x)
{
	return justi::Rg(AFAK(x));
}

double DieselMixture::k(double x, double T)
{
	return justi::k_Justi(T, AFAK(x));
}

bool DieselMixture::_initPrepertyTable()
{
	property = ArrayTable(8, 0);
	std::vector<std::string> arr =
	{ "已燃百分比",
		"缸内气体总质量","空气质量","废气质量","已喷燃油量",
		"广义过量空气系数","残余废气系数",
		"气体常数",
	};
	std::vector<std::string> arr1 =
	{
		"Fraction",
		"kg","kg","kg","kg",
		"","",
		"J/(kg·k)",
	};
	property.setParaName(arr);
	property.setParaUnit(arr1);

	for (double i = 0; i <= 1; i += 0.01)
	{
		std::vector<double> temp =
		{
			i,
			M_total(i),M_air(i),M_exhaust(i),gf * i,
			AFAK(i),M_exhaust(i) / M_total(i),
			Rg(i)
		};
		property.append(temp);
	}
	return true;
}
