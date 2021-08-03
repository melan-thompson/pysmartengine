#include "GasProperty.h"
double justi::Rg(double AFAK)//J.keenan & J.Kaye formula, used for ideal gas
{
	return 9.81 * (29.2647 - 0.0402 / AFAK);//J/(kg*K)
}

double justi::u_Justi(double T, double AFAK, double Tref)
{
	double temp1 = 489.6 + 46.4 / pow(AFAK, 0.93);
	double temp2 = 7.768 + 3.36 / pow(AFAK, 0.8);
	double temp3 = 0.0975 + 0.0485 / pow(AFAK, 0.75);
	double temp4 = 1356.8 + temp1 * (T - Tref) * 1.e-2 + temp2 * pow(T - Tref, 2) * 1.e-4 - temp3 * pow(T - Tref, 3) * 1.e-6;
	return 0.1445 * temp4 * 1.e3;//J/kg
}

double justi::h_Justi(double T, double AFAK, double Tref)
{
	return u_Justi(T, AFAK, Tref) + Rg(AFAK) * T;
}

double justi::cv_Justi(double T, double AFAK, double Tref)
{
	double h = 1.e-3;
	return (u_Justi(T + h, AFAK, Tref) - u_Justi(T - h, AFAK, Tref)) / 2 / h;//J/(kg*K)
}

double justi::M(double AFAK)//J.keenan & J.Kaye formula, used for ideal gas
{
	return 28.9705 + 0.0403 / AFAK;//g/mol
}



double justi::k_Justi(double T, double AFAK)
{
	return 1 + Rg(AFAK) / cv_Justi(T, AFAK);
}

double justi::cp_Justi(double T, double AFAK, double Tref)
{
	return cv_Justi(T, AFAK) * k_Justi(T, AFAK);
}

double justi::cv_mean_Justi(double T_begin, double T_end, double AFAK, double step)
{
	if (T_begin > T_end)
	{
		std::cout << "T_begin must be less than T_end when calculating mean cv!!!\n";
		exit(EXIT_FAILURE);
	}
	if (step <= 0)
	{
		std::cout << "step must be great than 0 when calculating mean cv!!!\n";
		exit(EXIT_FAILURE);
	}
	double sum = 0;
	for (double T = T_begin; T < T_end; T += step)
	{
		double cv1 = cv_Justi(T, AFAK); double cv2 = cv_Justi(T + step, AFAK);
		sum += (cv1 + cv2) / 2 * step;
	}
	return sum / (T_end - T_begin);
}

double justi::cp_mean_Justi(double T_begin, double T_end, double AFAK, double step)
{
	if (T_begin > T_end)
	{
		std::cout << "T_begin must be less than T_end when calculating mean cv!!!\n";
		exit(EXIT_FAILURE);
	}
	if (step <= 0)
	{
		std::cout << "step must be great than 0 when calculating mean cv!!!\n";
		exit(EXIT_FAILURE);
	}
	double sum = 0;
	for (double T = T_begin; T < T_end; T += step)
	{
		double cv1 = cp_Justi(T, 1.e8); double cv2 = cp_Justi(T + step, 1.e8);
		sum += (cv1 + cv2) / 2 * step;
	}
	return sum / (T_end - T_begin);
}

double justi::k_mean_Justi(double T_begin, double T_end, double AFAK, double step)
{
	return cp_mean_Justi(T_begin, T_end, AFAK, step) / cv_mean_Justi(T_begin, T_end, AFAK, step);
}

double schmidt::Mcv_pure_air(double T)//T(k),M*cv(kJ/(kmol*K))
{
	double a0 = 19.596;
	double b0 = 2.8773e-3;
	double c0 = -2.5407e-7;
	return a0 + b0 * T + c0 * T * T;//kJ/(kmol*K)
}

double schmidt::Mcv_pure_exhaust_gas(double T)//T(k),M*cv(kJ/(kmol*K))
{
	double ar = 19.893;
	double br = 5.020e-3;
	double cr = -5.959e-7;
	return ar + br * T + cr * T * T;
}

double FZacharias::Z(double p, double T, double AFAK)//采用F.Zacharias公式计算压缩因子，已验证其正确性
		//p(Pa),T(k),AFAK(1),Z(1)
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

double FZacharias::Rg(double p, double T, double AFAK)//used for high pressure, F.Zacharias formula
	//p(Pa),T(K),AFAK(1),Rg(J/(kg*K))
{
	double K = (AFAK - 1) / (0.0698 + AFAK);
	double A = 0.35 - 0.05 * pow(K, 0.765);
	double B = 11.1 + 14.3 * pow(K, 0.51);
	double C = 0.252 + 0.102 * pow(K, 0.401);
	return 9.81 * (29 + A + 0.1e-5 * p * (B / pow(T, C) - A));//J/(kg*K)
}
