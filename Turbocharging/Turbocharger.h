#ifndef TURBOCHARGER_H
#define TURBOCHARGER_H
#include "ParameterTable.h"
#include "ArrayTable.h"
#include "GasProperty\GasProperty.h"
#include <algorithm>



class CompressorMap
{
	ArrayTable MapData;
	void plot()
	{}
};

class Compressor
{
public:
	double n, pi, w, eff;
	ParameterTable CompressoreMap;
	ArrayTable record;

	Compressor() {}

	double massFlowRate(double P0, double PK, double T0, double WK, double eff)
	{
		/*double k = _k(_TK_, 1.e8);
		double Rg = _Rg(1.e8);
		double temp = (k - 1) / k;
		_TK_ = T0 * (1 + (pow(PK / P0, temp) - 1) / eff);
		return WK / (cp_Justi(_TK_, 1.e8) * (_TK_ - T0));*///((1. / temp) * Rg * T0 * (pow(PK / P0, temp) - 1) / eff);

		//先迭代求解等熵压缩后的温度
		//方程Tk=T0*(pi)^(k-1/k),k是(Tk,T0)的函数
		double Tk_iso = _TK_;
		auto fun = [=](double Tk)->double {double ktemp = justi::k_mean_Justi(T0, Tk); return Tk - T0 * pow(PK / P0, (ktemp - 1) / ktemp); };
		while (abs(fun(Tk_iso))>1.e-5)
		{
			double h = 0.01;
			double temp = (fun(Tk_iso + h) + fun(Tk_iso - h)) / (2 * h);
			Tk_iso = Tk_iso - fun(Tk_iso) / temp;
		}

		//再求解质量流量
		double mass_flow = WK / (eff * (justi::h_Justi(Tk_iso, 1.e8) - justi::h_Justi(T0, 1.e8)));

		//再求解实际的温度
		//方程WK=m*(h1-h0),h是温度的函数
		auto fun2 = [=](double Tk)->double {return WK - mass_flow * (justi::h_Justi(Tk, 1.e8) - justi::h_Justi(T0, 1.e8)); };
		while (abs(fun2(_TK_))>1.e-5)
		{
			double h = 0.01;
			double temp = (fun2(_TK_ + h) + fun2(_TK_ - h)) / (2 * h);
			_TK_ = _TK_ - fun2(_TK_) / temp;
		}
		return mass_flow;
	}

	double massFlowRate(double n, double pi)
	{

	}

	double efficiency(double n, double pi)
	{

	}
private:
	double _TK_;
};

class Turbo
{
	double n, pi, w, eff;
	ParameterTable TurboMap;
	Turbo(double Area) :_Area_(Area)
	{
		_flag_ = 0;
	}

	double massFlowRate(double N, double Pin, double Tin, double Pout, double Rin, double kin)
	{
		double piT = Pin / Pout;
		if (piT < 1)
		{
			return 0;
		}
		//有流量时，涡轮后的温度不确定，也用不到
		_DMT_ = _FLOW(N, Pin, Tin, Pout, 273.15, Rin, kin, 287, 1.4, _uTurbo(piT) * _Area_);
		_flag_ = 1;//流量计算完成
		return _DMT_;
	}

	double Work(double Tin, double kin, double Rin, double piT, double eff)//必须先算流量才能算功
	{
		if (_flag_ == 0)
		{
			std::cout << "Turbo calculation error,please update the mass flow first before you calculate the power!!!\n";
			system("pause");
			return 0;
		}
		double temp = kin / (kin - 1);
		_flag_ = 0;
		return eff * _DMT_ * Rin * Tin * temp * (1 - pow(piT, temp));
	}

private:
	double _Area_;
	double _DMT_;
	int _flag_;
	double _uTurbo(double piT)
	{
		return -0.443 + 1.8313 * piT - 0.7963 * pow(piT, 2) + 0.1508 * pow(piT, 3) - 0.01056 * pow(piT, 4);
	}

	double _FLOW(double N, double P1, double T1, double P2, double T2, double R1, double K1, double R2, double K2, double MUF)
	{
		if (MUF <= 1.e-12 || T1 <= 1.e-8 || T2 <= 1.e-8)
		{
			return 0;
		}

		double TMP2 = P2 / P1;

		if (TMP2 < 1.0)//从1流向2
		{
			double TMP1 = MUF / 6. / N * sqrt(2 * K1 / R1 / T1);

			double TMP3 = pow(2 / (K1 + 1), 1 / (K1 - 1));

			if (TMP2 > pow(TMP3, K1))
			{
				return TMP1 * P1 * sqrt((pow(TMP2, 2 / K1) - pow(TMP2, (K1 + 1) / K1)) / (K1 - 1));
			}
			else
			{
				return TMP1 * P1 / sqrt(K1 + 1) * TMP3;
			}
		}
		else
		{
			double TMP1 = MUF / 6 / N * sqrt(2 * K2 / R2 / T2);
			double TMP2 = P1 / P2;
			double TMP3 = pow(2 / (K2 + 1), 1 / (K2 - 1));
			if (TMP2 > pow(TMP3, K2))
			{
				return -TMP1 * P2 * sqrt((pow(TMP2, 2 / K2) - pow(TMP2, (K2 + 1) / K2)) / (K2 - 1));
			}
			else
			{
				return -TMP1 * P2 / sqrt(K2 + 1) * TMP3;
			}
		}
	}

};



#endif // !TURBOCHARGER_H

