#ifndef PIPE_H
#define PIPE_H
#include<cmath>
#define PI 3.141592654
#include "ParameterTable.h"
#include "ArrayTable.h"
#include "GasProperty.h"
#include <Eigen/Dense>


#ifndef VOLUME
#define VOLUME
class Volume//气缸的热物性
{
public:
	double V, T, M, AFA;
	double p;//pressure
	double Rg, k, cp, cv;//gas property, all by calculating
	double h, u;//焓和热力学能
	Volume() {};
	Volume(double Pressure, double Temperature, double Volume, double AFA) :p(Pressure), T(Temperature), V(Volume), AFA(AFA)
	{
		M = p * Volume / Rgfun(AFA) / T;
	};
	bool update()
	{
		k = _k(T, AFA);
		cv = cv_Justi(T, AFA);
		Rg = _Rg(AFA);
		cp = cp_Justi(T, AFA);
		u = u_Justi(T, AFA);
		p = M * Rg * T / V;
		while (true)
		{
			double p1 = _Z(p, T, AFA) * M * Rg * T / V;
			if (abs(p1 - p) / p < 1.e-5)
			{
				break;
			}
			p = p1;
		}
		h = u + p * V / M;
		return true;
	}

	double Rgfun(double AFA)
	{
		return _Rg(AFA);
	}

private:
	double _Z(double p, double T, double AFAK)
	{
		double T0 = T / 1000;
		double P0 = p * 1.e-6 / 0.098;
		double gamma = (AFAK - 1) / (1.0698 + AFAK);
		double A = 2.77105e-4 - 0.900711e-4 * gamma;
		double B = 6.42217e-4 - 0.98367e-4 * gamma;
		double D = 0.8868e-2 - 0.6131e-2 * gamma;
		double temp = -0.5 * gamma * pow(P0, -0.103) - 0.4 + (0.12 - 0.29 * pow(P0, -0.127)) * T0 / 1.65;
		double y1 = 15 * pow(T0 / 1.65, temp) - 6.9078 - ((1 - gamma) - 0.2772 * pow(1 - gamma, 4)) * 0.9088 * pow(P0, -0.0212) + (0.573 + 0.083 * log(P0)) * exp(-0.306 * T);
		double y2 = 18.0972 - 2.43 * T0 * pow(P0, -0.045);
		double temp1 = 1 + P0 / T0 * (B - A / T0 * exp(D / pow(T0, 2)));
		double temp2 = (0.420 - 0.193 * gamma) * pow(P0, 0.007) / (1 + 1 / T0 * exp(y1));
		double temp3 = (0.774 - 0.119 * gamma - (0.0128 + 0.005 * gamma) * log(P0)) / (1 + 1 / T0 * exp(y2));
		return temp1 + temp2 + temp3;
	}

	double _Rg(double AFAK)
	{
		return 9.81 * (29.2647 - 0.0402 / AFAK);//292.647 - 0.0402 / AFAK;//J/(kg*K)
	}

	double cv_Justi(double T, double AFAK = 10000, double Tref = 273.15)
	{
		double h = 1.e-3;
		return (u_Justi(T + h, AFAK, Tref) - u_Justi(T - h, AFAK, Tref)) / 2 / h;//J/(kg*K)
	}


	double _M(double AFAK)
	{
		return 28.9705 + 0.0403 / AFAK;//g/mol
	}

	double kploy(double T, double AFAK)
	{
		return 1.4373 - 1.318e-4 * T + 3.12 * 1.e-8 * T * T - 4.8e-2 / AFAK;
	}

	double RgZacharias(double p, double T, double AFAK)
	{
		double K = (AFAK - 1) / (0.0698 + AFAK);
		double A = 0.35 - 0.05 * pow(K, 0.765);
		double B = 11.1 + 14.3 * pow(K, 0.51);
		double C = 0.252 + 0.102 * pow(K, 0.401);
		return 9.81 * (29 + A + 0.1e-5 * p * (B / pow(T, C) - A));
	}

	double _k(double T, double AFAK)
	{
		return 1 + _Rg(AFAK) / cv_Justi(T, AFAK);
	}

	double cp_Justi(double T, double AFAK)
	{
		return cv_Justi(T, AFAK) * _k(T, AFAK);
	}

	double u_Justi(double T, double AFAK, double Tref = 273.15)
	{
		double temp1 = 489.6 + 46.4 / pow(AFAK, 0.93);
		double temp2 = 7.768 + 3.36 / pow(AFAK, 0.8);
		double temp3 = 0.0975 + 0.0485 / pow(AFAK, 0.75);
		double temp4 = 1356.8 + temp1 * (T - Tref) * 1.e-2 + temp2 * pow(T - Tref, 2) * 1.e-4 - temp3 * pow(T - Tref, 3) * 1.e-6;
		return 0.1445 * temp4 * 1.e3;//J/kg
	}
};

class Pipe :public Volume
{
	double M_exhaust()
	{
		return M * (1 + L0) / (L0 * AFA + 1);
	}

	double M_air()
	{
		return M - M_exhaust();
	}

	double M_burned()
	{
		return M_exhaust() / (1 + L0);
	}


	ArrayTable record;
	Pipe(int Type,double Vol,double pInit,double Tinit,double AFAInit,double TW,double heatexchangeCoeff=0):Type(Type),TW(TW),heatexchangeCoeff(heatexchangeCoeff)
	{
		V = Vol;
		p = pInit; T = Tinit; AFA = AFAInit;
		M = p * V / Rgfun(AFA) / T;
		L0 = 14.3;

		record = ArrayTable(16, 0);
		std::string paraName[] = { "曲轴转角","容积内压力","容积内温度","广义过量空气系数",
		"容积内质量","容积内废气质量","容积内空气质量","进入质量流量","出质量流量",
		"定容比热","定压比热","气体常数","比热容比",
		"传热率","进气焓流","排气焓流" };
		std::string paraUnit[] = { "CA","Pa","K","/",
			"kg","kg","kg","kg/s","kg/s",
			"J/(kg*K)","J/(kg*K)" ,"J/(kg*K)" ,"/",
			"W","W","W" };
		record.setParaName(paraName);
		record.setParaUnit(paraUnit);

	};

	bool RecordThisStep(double Fi)
	{
		double arr[] = { Fi,p,T,AFA,
			M,M_exhaust(),M_air(),_DMI_,_DMO_,
			cv,cp,Rg,k,
			_DQW_,_DHI_,_DHO_};
		record.append(arr);
		return true;
	}

	bool ClearRecord()
	{
		record.Clear();
	}


	Eigen::VectorXd D(double N,double S,double D,double DMI, double DMO, double DHI, double DHO, double AFAI, double AFAO)
	{
		_DHI_ = DHI; _DHO_ = DHO; _DMI_ = DMI; _DMO_ = DMO;
		update();//更新热力参数
		_DQW_ = _DQW(Type, S, D, N, TW);
		Eigen::VectorXd DXDFI(3);
		DXDFI(0) = _DM(DMI, DMO);
		DXDFI(1) = _DAFA(DMI, AFAI, DMO, AFAO);
		DXDFI(2) = _DT(DHI, DHO, _DQW_, DXDFI(0));
	}
private:
	int Type;
	double L0;
	double heatexchangeCoeff;
	double TW;
	double _DHI_, _DHO_, _DMI_, _DMO_, _DQW_;
	double _DM(double DMI, double DMO)
	{
		return DMI - DMO;
	}

	double _DT(double DHI, double DHO, double DQW, double DM)
	{
		return (DHI - DHO - DQW - u * DM) / cv / M;
	}

	double _DAFA(double DMI, double AFAI, double DMO, double AFAO)
	{
		double AFAIN, AFAOU;
		if (DMI > 0) AFAIN = AFAI;
		else
		{
			AFAIN = AFA;
		}
		if (DMO > 0)
		{
			AFAOU = AFA;
		}
		else
		{
			AFAOU = AFAO;
		}
		double TMP1 = (AFA * L0 + 1) / (AFAIN * L0 + 1);
		double TMP2 = (AFA * L0 + 1) / (AFAOU * L0 + 1);
		return (DMI * (1 - TMP1) - DMO * (1 - TMP2)) / M_burned() / L0;
	}

	double _DQW(int type,double S,double D,double N,double TW)
	{
		if (type==0)//无传热
		{
			return 0;
		}
		else if(type==1)//第一种传热方式
		{
			return heatexchangeCoeff * sqrt(p / 9.80665 * 1.e-4 * T) * pow(S * N / 30, 1. / 3) * pow(D, 2) * (0.5 + 0.78) * (T - TW) / 21600 / N;
		}
	}
};

#endif // !VOLUME
#endif // !PIPE_H

