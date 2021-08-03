#ifndef VALVE_H
#define VALVE_H
#include<cmath>
#define PI 3.141592654
#include "ParameterTable.h"
#include "ArrayTable.h"
#include "GasProperty/DieselMixtureProperty.h"

//一个喷管
class ValveSimple
{
public:
	double flowArea;
	ArrayTable flow;

	ValveSimple(double _flowArea = 1.);

	ArrayTable airFlowExample(double P1, double T1, double P2 = NULL, double T2 = NULL);

	double massFlowRate(double P1, double T1, double R1, double K1, double P2, double T2 = NULL, double R2 = NULL, double K2 = NULL);

private:

};

class Valve
{
public:

	double valveDiameter;//气阀
	double Open_Timing_Angle;
	double Close_Timing_Angle;
	double rock;//
	double sigma;
	ArrayTable Valve_Lift;
	ArrayTable record;
	ArrayTable flowCoeff;


	Valve() {};

	Valve(ArrayTable Valve_LiftIn, double valveDiameter, double Open_Timing_Angle, double Close_Timing_Angle, double SigmaIn, double rock = 1, ArrayTable _flow_coeff = NULL);
	
	~Valve() {};

	bool RecordThisStep(double Fi)
	{
		std::vector<double> arr = { Fi,_Lift_,_massFlow_ };
		record.append(arr);
		return true;
	}

	bool ClearRecord()
	{
		record.Clear();
	}

	double DM(double Fi, double N, double P1, double T1, double P2, double T2, double R1, double K1, double R2, double K2, double FlowCoeff = 1)
	{
		_massFlow_ = Flow(N, P1, T1, P2, T2, R1, K1, R2, K2, FlowCoeff * Flow_Area(Fi));
		return _massFlow_;
	};

	double Lift(double Fi)
	{
		_Lift_ = Valve_Lift.linearInterpolate(Fi, 1, 1);
		return _Lift_;
	}

	//气阀流量系数的经验公式
	double flowcoefficient(double _valveLift);
private:

	double _Lift_,_massFlow_;
	double Flow_Area(double Fi)
	{
		double HV = Valve_Lift.linearInterpolate(Fi, 1, 1);
		double U = 0.98 - 3.3 * pow((HV / valveDiameter), 2);
		return U * PI * HV * cos(sigma) * (valveDiameter + HV * sin(sigma) * cos(sigma));
	};//Valve Flow Area

	double Flow(double N, double P1, double T1, double P2, double T2, double R1, double K1, double R2, double K2, double MUF)
	{
		if (MUF <= 1.e-12 || T1 <= 1.e-8 || T2 <= 1.e-8)
		{
			return 0;
		}

		double TMP2 = P2 / P1;

		if (TMP2 < 1.0)
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
#endif // !VALVE_H


