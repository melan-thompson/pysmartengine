#include "Cylinder.h"
#define PI 3.141592654
#include "GasProperty.h"
#include "ParameterTable.h"
#include <assert.h>
#include <vector>

// ====================================================================================================
// Outflow routine for open end, partially open end, valve, or port; non-homentropic flow
//本子函数参照帝国理工论文第四章编写
// ====================================================================================================
double ondas_out(double lambda_in, double AA, double phi, double rb, double k, double tor)
{
	//初始特征值无量纲化
	double lambda_in_star = lambda_in / (AA * pow(rb, (k - 1) / (2 * k)));
	double lambda_out;
	// A_star lies between 1.0 and lambda_in_star - start at halfway
	double A_star = (lambda_in_star + 1) / 2;
	// Use a quarter of the difference the first time
	double del_A_star = (lambda_in_star - 1) / 4;
	double SUM;
	double A_star_prev;
	// Always halve step size even if there is no sign change in SUM (leads to less loops)
	bool ALWAYS_HALVE = true;
	bool GREATER = false;
	int loop = 0;
	//迭代计算音速，当音速满足收敛条件时退出迭代
	do {
		//记录初始音速
		A_star_prev = A_star;
		++loop;
		// Eqn (7.122)
		SUM = (pow(A_star, 4 / (k - 1)) - pow(phi, 2)) * pow((lambda_in_star - A_star), 2) - ((k - 1) / 2) * pow(phi, 2) * (pow(A_star, 2) - 1);
		if (SUM < 0.0) {
			A_star -= del_A_star;
			if (ALWAYS_HALVE || (GREATER && loop > 1)) del_A_star /= 2;
			GREATER = false;
		}
		else {
			A_star += del_A_star;
			if (ALWAYS_HALVE || (!GREATER && loop > 1)) del_A_star /= 2;
			GREATER = true;
		}
		//}while(fabs(SUM)>pPpt->NOZZ_TOL);
	} while (fabs(SUM) > tor || fabs(A_star - A_star_prev) > tor);

	// Using eqn (7.118)
	double lambda_out_star = 2 * A_star - lambda_in_star;
	// Using eqn (7.117)
	double U_star_temp = (lambda_in_star - lambda_out_star) / (k - 1);
	// Eqn (7.124)
	double Mt = U_star_temp * pow(A_star, 2 / (k - 1)) / phi;

	if (Mt >= 1.0/*choked flow*/ || Mt < 0/*unrealistic flow situation: inflow during outflow routine*/) {
		// From eqn (7.128) upper limit of A_star_cr is ((k+1)/2)^0.5, lower limit is 1, start at halfway
		double A_star_cr = (1 + sqrt((k + 1) / 2)) / 2;
		// Start with 1/4 of the difference
		double del_A_star_cr = (sqrt((k + 1) / 2) - 1) / 4;
		double SUM_cr;
		GREATER = false;
		loop = 0;
		do {
			++loop;
			// Eqn (7.128)
			SUM_cr = pow(phi, 2) - ((k + 1) / (k - 1) - (2 / (k - 1)) * pow(A_star_cr, 2)) * pow(A_star_cr, 4 / (k - 1));
			if (SUM_cr > 0.0) {
				A_star_cr -= del_A_star_cr;
				if (ALWAYS_HALVE || (!GREATER && loop > 1)) del_A_star_cr /= 2;
				GREATER = true;
			}
			else {
				A_star_cr += del_A_star_cr;
				if (ALWAYS_HALVE || (GREATER && loop > 1)) del_A_star_cr /= 2;
				GREATER = false;
			}
		} while (fabs(SUM_cr) > tor);
		// Using eqn (7.131)
		lambda_out = (1 - ((k - 1) / 2) * phi * pow(1 / A_star_cr, (k + 1) / (k - 1))) / (1 + ((k - 1) / 2) * phi * pow(1 / A_star_cr, (k + 1) / (k - 1))) * lambda_in;
	}
	else {
		lambda_out = lambda_out_star * AA * pow(rb, (k - 1) / (2 * k));
	}

	//返回特征值
	return lambda_out;
}

// ====================================================================================================
// Inflow routine for open end, partially open end, valve, or port; non-homentropic flow
// ====================================================================================================
std::vector<double> ondas_valve(double clin, double clout, double AAn, double psi, double pp, double A0, double k, bool INFLOW_CONST_P, double tolLoop)
///
{

	double Ac, rc, AA_is, del_AA;
	double U, C, R, PpbyPc, AA_calc;
	double lambda_in_c, lambda_out_c, AA_est, lambda_in_c_prev, lambda_out_c_prev;
	bool CHOKED, SONIC;
	bool GREATER, GREATER_PREV;
	bool RESET;

	GREATER = false;
	GREATER_PREV = false;



	// Numbers in braces refer to steps in Benson's (1982, pg. 378) non-homentropic outflow to a pipe routine (Benson, 1982, pg. 378)
	// ----------------------------------------------------------------------------------------------------

	// (1) Enter with known characteristic values from boundary point
	// ----------------------------------------------------------------------------------------------------
	Ac = A0;	// Non-dimensional speed of sound in source (e.g., cylinder)

	// (2) Set uncorrected values
	// ----------------------------------------------------------------------------------------------------
	// Entropy level for isentropic process
	AA_is = Ac / pow(pp, (k - 1) / (2 * k));
	// Use uncorrected value as initial estimate
	lambda_in_c = clin;
	// Use previous result fot initial estimate
	lambda_out_c = clout;

	// (3) Test rc
	// ----------------------------------------------------------------------------------------------------
	rc = pp;
	if (rc < 1.0) { AA_est = AA_is; del_AA = 0.5; }
	else {
		AA_est = AA_is;
		del_AA = fabs(Ac - AA_is) / 2;
		if (del_AA == 0) del_AA = 1e-6;
	}


	int loop = 0;
	do {
		RESET = false;
		++loop;

		// (4) Calculate new value of lambda_in from eqn (7.155)
		// ----------------------------------------------------------------------------------------------------
		lambda_in_c_prev = lambda_in_c;
		//lambda_in_c = 2*lambda_in_n*(AA_est/(AA_est + AA_n)) + lambda_out_c*((AA_est - AA_n)/(AA_est + AA_n));
		lambda_in_c = ((lambda_in_c + lambda_out_c) / 2) * (1 - ((AAn) / AA_est)) + (clin); // Eqn 7.155

		// (5) Test (Ac - lambda_in_c)
		// ----------------------------------------------------------------------------------------------------
		if (Ac - lambda_in_c < 0) {
			RESET = true;
			AA_est = AA_is;
			del_AA /= 2;
			lambda_in_c_prev = lambda_in_c;
			//lambda_in_c = 2*lambda_in_n*(AA_est/(AA_est + AA_n)) + lambda_out_c*((AA_est - AA_n)/(AA_est + AA_n));
			lambda_in_c = ((lambda_in_c + lambda_out_c) / 2) * (1 - ((AAn) / AA_est)) + (clin); // Eqn 7.155
		}

		// (6) Calculate new value of lambda_out from eqn (7.154)
		// ----------------------------------------------------------------------------------------------------
		lambda_out_c_prev = lambda_out_c;
		lambda_out_c = ((3 - k) / (k + 1)) * lambda_in_c + (2 / (k + 1)) * sqrt((pow(k, 2) - 1) * pow(Ac, 2) + 2 * (1 - k) * pow(lambda_in_c, 2));

		// (7) Calculate U from eqn (7.149)
		// ----------------------------------------------------------------------------------------------------
		U = (lambda_out_c - lambda_in_c) / ((k - 1) * Ac);

		if (INFLOW_CONST_P) { // Inflow constant pressure (valve) model,等压模型

			// (8) Calculate C from eqn (7.139)
			// ----------------------------------------------------------------------------------------------------
			C = (((k - 1) / 2) * pow(U, 2)) / pow((1 - ((k - 1) / 2) * pow(U, 2)), 2);

			// (9) Test 4C/(k^2 - 1) - psi^2
			// ----------------------------------------------------------------------------------------------------
			if ((4 * C) / (pow(k, 2) - 1) - pow(psi, 2) >= 0) {
				// Sonic flow in the valve
//cout << "Sonic flow in the valve" << endl;
				// (10) If sonic flow, calculate pp/pc from eqn (7.142)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = true;
				SONIC = true;
				if (U == 0) PpbyPc = 1.0;
				else PpbyPc = psi * pow(2 / (k + 1), (k + 1) / (2 * (k - 1))) * (1 - ((k - 1) / 2) * pow(U, 2)) * (1 / U);
				// Can happen if U is very small but not zero
				if (PpbyPc == 0) PpbyPc = 1;
			}
			else {
				// (10) If subsonic flow calculate pp/pc fom eqn (7.138)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = false;
				SONIC = false;
				if (U == 0) PpbyPc = 1.0;
				else PpbyPc = pow((1 / (2 * C)) * (psi * sqrt(pow(psi, 2) + 4 * C) - pow(psi, 2)), k / (k - 1));
				if (PpbyPc == 0) PpbyPc = 1; // Can happen if U is very small but not zero
			}
		}
		else {
			// Inflow pressure drop (port) model
			// (8) Calculate R from eqn (7.157)
			// ----------------------------------------------------------------------------------------------------
			R = (pow(psi, 2) * pow(1 - ((k - 1) / 2) * pow(U, 2), (k - 3) / (k - 1))) / ((k - 1) * pow(U, 2));

			// (9) Test (2/(k+1))*pow(1 - ((k-1)/2)*pow(U,2), 2/(k-1))*pow(U,2) - pow(psi,2)
			// ----------------------------------------------------------------------------------------------------
			if ((2 / (k + 1)) * pow(1 - ((k - 1) / 2) * pow(U, 2), 2 / (k - 1)) * pow(U, 2) - pow(psi, 2) >= 0) {
				// Sonic flow in the ports
				// (10) If sonic flow, calculate pp/pc from eqn (7.159)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = true;
				SONIC = true;
				if (U == 0) PpbyPc = 1;
				else PpbyPc = psi * pow(2 / (k + 1), (k + 1) / (2 * (k - 1))) * (1 - ((k - 1) / 2) * pow(U, 2)) * (1 / U);
				if (PpbyPc == 0) PpbyPc = 1; // Can happen if U is very small but not zero
			}
			else {
				// (10) If subsonic flow calculate pp/pc fom eqn (7.156)
				// ----------------------------------------------------------------------------------------------------
				CHOKED = false;
				SONIC = false;
				if (U == 0) PpbyPc = 1;
				else PpbyPc = pow(R * (pow(1 + (2 / R) * (1 - ((k - 1) / 2) * pow(U, 2)), 0.5) - 1), k / (k - 1));
				if (PpbyPc == 0) PpbyPc = 1; // Can happen if U is very small but not zero
			}
		}

		// (11) Calculate AA from eqn (7.153)
		// ----------------------------------------------------------------------------------------------------
		AA_calc = ((lambda_in_c + lambda_out_c) / 2) * pow(1 / PpbyPc, (k - 1) / (2 * k)) * pow(1 / rc, (k - 1) / (2 * k));

		// (12) Test {(AA)calc - AAest}
		// ----------------------------------------------------------------------------------------------------
		if (loop > 1) GREATER_PREV = GREATER;
		if (AA_calc - AA_est > 0) {
			AA_est += del_AA;
			GREATER = true;
		}
		else {
			AA_est -= del_AA;
			GREATER = false;
		}
		if (GREATER_PREV != GREATER && loop > 1) del_AA /= 2;
	} while (fabs(AA_calc - AA_est) >= tolLoop || fabs(lambda_in_c - lambda_in_c_prev) >= tolLoop || fabs(lambda_out_c - lambda_out_c_prev) >= tolLoop);

	//返回计算结果，更新特征值
	/**clin = lambda_in_c;
	*clout = lambda_out_c;
	*AAn = AA_est;*/

	return std::vector<double> {lambda_in_c, lambda_out_c, AA_est};
}


heatRelease::SingleWiebe::SingleWiebe(double SOC, double TCD, double m, double a)
	:SOC(SOC), TCD(TCD), m(m), a(a)
{
	for (double i = SOC-10; i < SOC + TCD+10; i += 0.1)
	{
		std::vector<double> arr = { i,DX(i),X(i) };
		data.append(arr);
	}
}

double heatRelease::SingleWiebe::X(double CA)
{
	if (CA < SOC)return 0;
	double temp = -a * pow((CA - SOC) / TCD, m + 1);
	return 1 - exp(temp);
}

double heatRelease::SingleWiebe::DX(double CA)
{
	if (CA < SOC)return 0;
	double temp = -a * pow((CA - SOC) / TCD, m + 1);
	return a * (m + 1) / TCD * pow((CA - SOC) / TCD, m) * exp(temp);
}

double heatRelease::SingleWiebe::getFractionCA(double fraction)
{
	if (fraction < 0 || fraction>1) throw std::range_error("burn mass fraction should be in [0,1]");
	return data.linearInterpolate(fraction,0,2);
}


heatRelease::DoubleWiebe::DoubleWiebe(double PF, double PSOC, double PTCD, double DTCD, double PM, double DM, double DSOCin)
	:PF(PF), PSOC(PSOC), PTCD(PTCD), PM(PM), DSOC(DSOCin), DTCD(DTCD), DM(DM)
{
	if (PF < 0 || PF>1)
	{
		std::cout << "Premixed combustion fraction value should between 0 to 1";
		system("pause");
		exit(EXIT_FAILURE);
	}
	DF = 1 - PF;
	if (DSOCin == INFINITY)
	{
		DSOC = PSOC + PTCD / 2;
	}
	PremixWiebe = SingleWiebe(PSOC, PTCD, PM);
	DisfusionWiebe = SingleWiebe(DSOC, DTCD, DM);

	for (double i = PSOC; i < DSOC + DTCD; i += 0.1)
	{
		data.append(std::vector<double>{ i, DX(i), X(i) });
	}
}

double heatRelease::DoubleWiebe::X(double CA)
{
	return PremixWiebe.X(CA) * PF + DisfusionWiebe.X(CA) * DF;
}

double heatRelease::DoubleWiebe::DX(double CA)
{
	return PremixWiebe.DX(CA) * PF + DisfusionWiebe.DX(CA) * DF;
}

double heatRelease::DoubleWiebe::getFractionCA(double fraction)
{
	if (fraction < 0 || fraction>1) throw std::range_error("burn mass fraction should be in [0,1]");
	return data.linearInterpolate(fraction, 0, 2);
}


heatRelease::HeatReleaseData::HeatReleaseData(ArrayTable data)
{
	ArrayTable temp(3, 0);
	temp.table[0].ColName = "Crank angle";
	temp.table[0].ColUnit = "CA";
	temp.table[1].ColName = "DX";
	temp.table[2].ColName = "X";
	temp.table[1].ColUnit = "1/deg";
	for (int i = 0; i < data.row; i++)
	{
		std::vector<double> arr = { data.table[0].data[i],data.table[1].data[i],data.Integrate(1,0,i) };
		temp.append(arr);
	}

	std::cout << "Total burn fraction is " << temp.table[2].data[temp.row - 1];
	heat_release_data = temp;
};


double heatRelease::HeatReleaseData::SOC()//Start of combustion
{
	return heat_release_data.table[0].data[findSOCIndex()];
}

double heatRelease::HeatReleaseData::EOC()//End of combustion
{

	return heat_release_data.table[0].data[findEOCIndex()];

}

double heatRelease::HeatReleaseData::DX(double Fi)
{
	return heat_release_data.linearInterpolate(Fi, 1, 1);
}

double heatRelease::HeatReleaseData::X(double Fi)
{
	return heat_release_data.linearInterpolate(Fi, 2, 1);
}

ParameterTable heatRelease::HeatReleaseData::analyze()
{
	ParameterTable result;
	result.append("Start of combustion", "°CA", SOC());
	result.append("End of combustion", "°CA", EOC());
	double CA10 = NULL, CA50 = NULL, CA90 = NULL;
	for (int i = 0; i < heat_release_data.row; i++)
	{
		if (heat_release_data.table[2].data[i] > 0.1)
		{
			CA10 = heat_release_data.table[0].data[i];
			break;
		}
	}

	for (int i = 0; i < heat_release_data.row; i++)
	{
		if (heat_release_data.table[2].data[i] > 0.5)
		{
			CA50 = heat_release_data.table[0].data[i];
			break;
		}
	}

	for (int i = 0; i < heat_release_data.row; i++)
	{
		if (heat_release_data.table[2].data[i] > 0.9)
		{
			CA90 = heat_release_data.table[0].data[i];
			break;
		}
	}

	result.append("Start of combustion", "°CA", SOC());
	result.append("End of combustion", "°CA", EOC());
	result.append("Combustion duration", "°CA", EOC() - SOC());
	result.append("CA10", "°CA", CA10);
	result.append("CA50", "°CA", CA50);
	result.append("CA90", "°CA", CA90);
	return result;
}

int heatRelease::HeatReleaseData::findSOCIndex()//Start of combustion
{
	for (int i = 0; i < heat_release_data.row; i++)
	{
		if (heat_release_data.table[1].data[i] > 0.000001)
			return i;
	}
	std::cout << "Heat release data has some problem!!!Can not find the start of combustion point!!!Will return -inf";
	system("pause");
	return -INFINITY;
}

int heatRelease::HeatReleaseData::findEOCIndex()//End of combustion
{
	for (int i = heat_release_data.row - 1; i > 0; i--)
	{
		if (heat_release_data.table[1].data[i] > 0.000001)
			return i;
	}
	std::cout << "Heat release data has some problem!!!Can not find the end of combustion point!!! Will return -inf";
	system("pause");
	exit(EXIT_FAILURE);
	return -INFINITY;
}

double heatTransfer::CylinderHeatTransfer::HeatTrasferRate(double S, double D, double V, double N, double P, double T)
{
	double CM = N * S / 30;
	double AL = 4. * V / D;
	double AFW = CQW * pow(D, -0.214) * pow(CM * P / 9.80665 * 1.e-4, 0.786) * pow(T, -0.525);
	return AFW * (AL * (T - CT) + PI * pow(D, 2) / 4. * (2 * T - HT - PT)) / 21600 / N;
};

heatTransfer::CylinderHeatTransfer::CylinderHeatTransfer(double Head_Temperature, double Piston_Temperature, double Cylinder_Temperature, double CorrectCoeff) :HT(Head_Temperature), PT(Piston_Temperature), CT(Cylinder_Temperature)
{
	CQW = 1256040 * CorrectCoeff;
}

//Pe为平均有效压力
heatTransfer::CylinderHeatTransfer::CylinderHeatTransfer(double Pe, double CorrectCoeff)
{
	CQW = 1256040 * CorrectCoeff;
	HT = 100 + 70 * Pe / 1.e5;
	PT = 100 + 120 * Pe / 1.e5;
	CT = 100 + 40 * Pe / 1.e5;
}

///<summary>Sitkei传热公式,适用于小型柴油机计算</summary>
double heatTransfer::Sitkei(double p, double T, double Cm, double D, double S, double V)
{
	double h = V / (PI / 4 * D * D);//活塞顶面到气缸盖燃烧室表面的距离
	double de = (2 * D * h) / (D + 2 * h);
	double b = 0.25;
	//对于直喷式燃烧室b=0.0-0.15
	//涡流室式燃烧室b=0.15-0.30
	//预燃室式燃烧室b=0.25-0.40
	double ag = 1.294e-2 * (1 + b) * pow(de, -0.3) * pow(T, 0.2) * pow(p, 0.2) * pow(Cm, 0.7);//(W/(m^2*K))
	return ag;
}

///<summary>Woschni公式</summary>
double heatTransfer::simpleWoschini(double D, double vm, double pz, double Tz)
{
	return 0.303 * pow(D, -0.214) * pow(vm * pz / 1.e6, 0.786) * pow(Tz, -0.525);
}

CylinderGeometry::CylinderGeometry(double Bore, double Stroke, double CR, double length) :D(Bore), S(Stroke), Eps(CR)
{
	Lam = S / 2 / length;
	STRC2 = PI * pow(D, 2) / 4. * (S / (Eps - 1.) + S / 2. + S / Lam / 2.);
	STRC3 = PI * pow(D, 2) * S / 4. / 2.;
	STRC4 = STRC3 / Lam;
	STRC5 = STRC3 * PI / 180;
	VH = PI / 4 * pow(D, 2) * S;

	ArrayTable temp(2, 0);
	for (double i = 0; i < 720; i += 0.5)
	{
		std::vector<double> arr = { i,V(i) };
		temp.append(arr);
	}
	data = temp;
};

double CylinderGeometry::V(double Fi) const
{
	double Temp = Fi * PI / 180.;
	double H1 = cos(Temp), H2 = sin(Temp), H3 = sqrt(1 - pow(Lam, 2) * pow(H2, 2));
	return STRC2 - STRC3 * H1 - STRC4 * H3;
}

double CylinderGeometry::DV(double Fi)const
{
	double Temp = Fi * PI / 180;
	double H1 = cos(Temp), H2 = sin(Temp), H3 = sqrt(1 - pow(Lam, 2) * pow(H2, 2));
	return STRC5 * H2 * (1 + Lam * H1 / H3);
}

double CylinderGeometry::clearanceVolume()const
{
	return V(0);
}

double CylinderGeometry::TDCclearanceHeight()const
{
	return clearanceVolume() / (PI / 4 * pow(D, 2));
}

double CylinderGeometry::heatExchangeArea(double crank_angle)const
{
	double area = PI / 4 * pow(D, 2);
	return 2 * area + 4 * V(crank_angle) / D;
}

void CylinderGeometry::plotVolume()
{
	data.plot();
}

Injection::Injection(double MF, double L0, double QH) :MF(MF), L0(L0), QH(QH) {};


Volume::Volume(double Pressure, double Temperature, double Volume, double AFA) :p(Pressure), T(Temperature), V(Volume), AFA(AFA)
{
	if (AFA < 1)throw std::range_error("alphak can not be less than 1");
	M = p * Volume / _Rg(AFA) / T;
	data = ArrayTable(11, 0);
	updateStatus();

	/*std::vector<std::string> temp = {"Pressure","Temperature"};*/
	std::vector< std::string> paraName = {"time","pressure","temperature","alphak",
		"m","exhasut gas mass","air mass",
		"cp","cv","Rg","k",
		};//热力学能,//已燃百分比//
	std::vector< std::string> paraUnit = {"s","Pa","K","/",
		"kg","kg","kg",
		"J/(kg*K)","J/(kg*K)" ,"J/(kg*K)" ,"/", };
	std::vector<double> temp = {0,p,T,AFA,M,M_exhasut(),M_air(),cp,cv,Rg,k };
	data.setParaName(paraName);
	data.setParaUnit(paraUnit);
	data.append(temp);

};


bool Volume::updateStatus()
{
	k = _k(T, AFA); Rg = _Rg(AFA);
	cv = cv_Justi(T, AFA); cp = cp_Justi(T, AFA);
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
	u = u_Justi(T, AFA); h = u + p * V / M;
	return true;
}

void Volume::recordThisStep(double t)
{
	std::vector<double> temp = {t,p,T,AFA,M,M_exhasut(),M_air(),cp,cv,Rg,k };
	data.append(temp);
}

void update()
{
	
}

void Volume::showProp()
{
	data.show();
}

double Volume::DMDT()
{
	double result = 0;
	//加上上游流入的流量
	for (auto i = flowLast.begin(); i != flowLast.end(); ++i)
	{
		//std::cout << "Upflow" << (*i)->DMDT;
		result += (*i)->DMDT;
	};

	//减去下游流入的流量
	for (auto i = flowNext.begin(); i !=flowNext.end(); ++i)
	{
		//std::cout << "down flow" << (*i)->DMDT;
		result -= (*i)->DMDT;
	}

	//将质量变化记录下来，以供DTDT使用
	dm = result;
	return result;
}

double Volume::DTDT(double DQDT)
{
	double Hin = 0.f;
	for (auto i = flowLast.begin(); i != flowLast.end(); ++i)//遍历所有上游节点
	{
		if ((*i)->DMDT >=0)//从上游流入,为上游节点的焓乘以流量
		{
			//计算上游节点流入的焓
			//是否是调用的volume类的成员函数？
			Hin += (*i)->comLast->get_h()* (*i)->DMDT;
		}
		else//流出，为本节点的焓乘以流量
		{
			Hin += h * (*i)->DMDT;
		}
	}

	double Hout = 0.f;
	for (auto i = flowNext.begin(); i != flowNext.end(); ++i)
	{
		if ((*i)->DMDT >0)//从本节点流出，焓减少
		{
			Hout += h * (*i)->DMDT;
		}
		else//从下游节点流入，焓增加
		{
			Hout += (*i)->comLast->get_h() * (*i)->DMDT;
		}

	}

	//忽略成分的变化来计算温度
	double result = (Hin-Hout-DQDT-u*dm)/M/cv;

	return result;
}

double Volume::DAFA()
{
	double result = 0.f;
	for (auto i = flowLast.begin(); i != flowLast.end(); ++i)
	{
		if ((*i)->DMDT >0)
		{
			result += (*i)->DMDT * (1.f - (L0 * AFA + 1.f) / (L0 * (*i)->comLast->get_AFA()+1.f));
		}
	}

	for (auto i = flowNext.begin(); i != flowNext.end(); ++i)
	{
		if ((*i)->DMDT <0)//如果反向流动则容积内的成分会发生变化，否则不会，如果dmout<0,且下游更浓，则加浓
		{
			result -= (*i)->DMDT * (1.f - (L0 * AFA + 1.f) / (L0 * (*i)->comNext->get_AFA() + 1.f));
		}

	}

	//return result/(L0*M_exhasut()/(1+L0));

	return result / (L0 * M / (1 + L0 + AFA));
}

double Volume::M_exhasut()
{
	return M*(1+L0)/(L0*AFA+1);
}

double Volume::M_air()
{
	return M*L0*(AFA-1)/(L0*AFA+1);
}


double Volume::_Z(double p, double T, double AFAK)//已验证其正确性
{
	return FZacharias::Z(p, T, AFAK);
}

double Volume::_Rg(double AFAK)
{
	return justi::Rg(AFAK);
}

double Volume::cv_Justi(double T, double AFAK, double Tref)
{
	return justi::cv_Justi(T, AFAK, Tref);
}

double Volume::_M(double AFAK)
{
	return justi::M(AFAK);
}


double Volume::_k(double T, double AFAK)
{
	return justi::k_Justi(T, AFAK);
}

double Volume::cp_Justi(double T, double AFAK)
{
	return justi::cp_Justi(T, AFAK);
}

double Volume::u_Justi(double T, double AFAK, double Tref)
{
	return justi::u_Justi(T, AFAK, Tref);
}

Cylinder::Cylinder(double Fi, double pin, double Tin, double AFAin, double N, CylinderGeometry Geo, heatRelease::HeatReleaseData HRR, heatTransfer::CylinderHeatTransfer CylHT, Injection Inj) :Fi(Fi), Geo(Geo), HRR(HRR), N(N), CylHT(CylHT), Inj(Inj)//初始化缸内参数，曲轴转角，压力，温度
{
	//using namespace gasproperty;
	p = pin; T = Tin; AFA = AFAin;
	M = p * Geo.V(Fi) / justi::Rg(AFA) / T;
	Xk = M / (1 + Inj.L0 * AFA) / Inj.MF;
	V = Geo.V(Fi);
	DX = HRR.DX(Fi);

	record = ArrayTable(19, 0);
	std::vector< std::string> paraName = { "曲轴转角","气缸容积","缸内压力","缸内温度","广义过量空气系数",
		"缸内质量","废气质量","空气质量","进入质量流量","出质量流量",
		"定容比热","定压比热","气体常数","比热容比",
		"燃烧放热率","气缸传热率","进气焓流","排气焓流","指示功" };//热力学能,//已燃百分比//
	std::vector< std::string> paraUnit = { "CA","m^3","Pa","K","/",
		"kg","kg","kg","kg/s","kg/s",
		"J/(kg*K)","J/(kg*K)" ,"J/(kg*K)" ,"/",
		"J/deg","J/deg","J/deg","J/deg","J/deg" };
	record.setParaName(paraName);
	record.setParaUnit(paraUnit);
}

bool Cylinder::RecordThisStep()
{
	std::vector<double> arr = { Fi, V,p,T,AFA,
		M,M_exhaust(),M_air(),_DMI_,_DMO_,
		cv,cp,Rg,k,
		_DQ_,_DQW_,_DHI_,_DHO_,_DW_ };
	record.append(arr);
	return true;
}

bool Cylinder::ClearRecord()
{
	record.Clear();
	return true;
}


std::vector<double> Cylinder::D(double DMI, double DMO, double DHI, double DHO, double AFAI, double AFAO)
{
	_DMI_ = DMI; _DMO_ = DMO; _DHI_ = DHI; _DHO_ = DHO;
	V = Geo.V(Fi); updateStatus();//更新缸内的参数
	///////T,AFA,M,Fi由外部更新,V,等其他参数由内部更新/////

	_DW_ = p * Geo.DV(Fi);
	Xk = M / (1 + Inj.L0 * AFA) / Inj.MF;
	DX = HRR.DX(Fi);
	_DQ_ = DX * Inj.MF * Inj.QH;
	_DQW_ = _DQW();

	std::vector<double> DXDFI(3);
	DXDFI[0] = _DM(DMI, DMO);
	DXDFI[1] = _DAFA(DMI, AFAI, DMO, AFAO, DX);
	DXDFI[2] = _DT(DHI, DHO, _DQW_, _DQ_, _DW_, DXDFI[0]);//DT(DHI, DHO, DXDFI(0));
	return DXDFI;
};


double Cylinder::M_exhaust()
{
	return M * (1 + Inj.L0) / (Inj.L0 * AFA + 1);
}

double Cylinder::M_air()
{
	return M - M_exhaust();
}


double Cylinder::_DM(double DMI, double DMO)
{
	return DMI - DMO + Inj.MF * DX;
};

double Cylinder::_DT(double DHI, double DHO, double DQW, double DQ, double DW, double DM)
{
	return (DQ + DHI - DHO - DQW - DW - u * DM) / cv / M;
};

double Cylinder::_DMSC(double DMI, double DMO, double AFAOUT)
{
	if (abs(DMI) > 0)
	{
		if (DMO > 0)
		{
			return DMI - DMO * (AFA - 1.) * Inj.L0 / (AFA * Inj.L0 + 1);
		}
		else
		{
			return DMI - DMO * (AFAOUT - 1.) * Inj.L0 / (AFAOUT * Inj.L0 + 1);
		}
	}
	else
	{
		return 0;
	}

}

double Cylinder::_DAFA(double DMI, double AFAI, double DMO, double AFAO, double DX)
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
	double TMP1 = (AFA * Inj.L0 + 1) / (AFAIN * Inj.L0 + 1);
	double TMP2 = (AFA * Inj.L0 + 1) / (AFAOU * Inj.L0 + 1);
	return (DMI * (1 - TMP1) - DMO * (1 - TMP2) - DX * (Inj.L0 * AFA + 1) * Inj.MF) / Inj.MF / Xk / Inj.L0;
}

double Cylinder::_DQW()
{
	return CylHT.HeatTrasferRate(Geo.S, Geo.D, V, N, p, T);
}


bool CylinderPressure::readInFile(std::string filename)
{
	data.readInFile(filename);
	return true;
}

ArrayTable CylinderPressure::doMultipointSmooth(smoothtype type)
{
	ArrayTable AfterSmoothTable(2, data.row);
	if (type == smoothtype::_5point3times)
	{
		std::vector<double> a = data.table[1].data;
		std::vector<double> result(data.row);

		int i = -1;
		result[i + 1] = (69. * a[i + 1] + 4. * (a[i + 2] + a[i + 4]) - 6. * a[i + 3] - a[i + 5]) / 70.;
		result[i + 2] = (2. * (a[i + 1] + a[i + 5]) + 27. * a[i + 2] + 12. * a[i + 3] - 8 * a[i + 4]) / 35.;
		for (int i = 2; i < data.row - 2; i++)
		{
			result[i] = (-3. * (a[i - 2] + a[i + 2]) + 12. * (a[i - 1] + a[i + 1]) + 17. * a[i]) / 35.;
		}
		i = data.row - 1;
		result[i - 1] = (2. * (a[i - 4] + a[i]) - 8. * a[i - 3] + 12. * a[i - 2] - 27 * a[i - 1]) / 35.;
		result[i] = (-a[i - 4] + 4. * (a[i - 3] + a[i - 1]) - 6. * a[i - 2] + 69. * a[i]) / 70;

		AfterSmoothTable.table[0] = data.table[0];
		AfterSmoothTable.table[1].ColName = data.table[1].ColName;
		AfterSmoothTable.table[1].ColUnit = data.table[1].ColUnit;
		AfterSmoothTable.table[1].data = result;
	}
	else
	{
		for (int i = 0; i < data.row; i++)
		{
			AfterSmoothTable.table[0].data[i] = data.table[0].data[i];
			AfterSmoothTable.table[1].data[i] = MultipointSmooth(i, type, data);
		}
	}
	dataAfterSmooth = AfterSmoothTable;
	return AfterSmoothTable;
}

double CylinderPressure::findStartOfCombustion()
{
	data.diff(1); data.diff(2); data.diff(3);
	int startcom = data.findMaxDataIndex(4);
	return data.table[0].data[startcom];
}

void CylinderPressure::plot()
{
	data.plot(1);
}

ArrayTable CylinderPressure::LogPLogV(CylinderGeometry CylGeo)
{
	ArrayTable temp(2, 0);
	for (int i = 0; i < data.row; i++)
	{
		std::vector<double> arr = { log(CylGeo.V(data.table[0].data[i])),log(data.table[1].data[i]) };
		temp.append(arr);
	}
	return temp;
}

ArrayTable CylinderPressure::PVDiagram(CylinderGeometry CylGeo)
{
	ArrayTable temp(2, 0);
	for (int i = 0; i < data.row; i++)
	{
		std::vector<double> arr = { CylGeo.V(data.table[0].data[i]),data.table[1].data[i] };
		temp.append(arr);
	}
	return temp;
}



ArrayTable CylinderPressure::netHeatReleaseRate(CylinderGeometry CylGeo, DieselMixture mix)
{
	auto gamma = [](double T)->double {return 1.338 - 6.0e-5 * T + 1.0e-8 * T * T; };
	ArrayTable temp(2, 0);
	ArrayTable datacopy = data;
	datacopy.diff(1);
	for (int i = 0; i < datacopy.row; i++)
	{
		double Ttemp = datacopy.table[1].data[i] * CylGeo.V(datacopy.table[0].data[i]) / mix.Rg(0) / mix.M_total(0);
		double k = gamma(Ttemp);
		double dQ = k / (k - 1) * datacopy.table[1].data[i] * CylGeo.DV(datacopy.table[0].data[i]) + 1 / (k - 1) * CylGeo.V(datacopy.table[0].data[i]) * datacopy.table[2].data[i];
		std::vector<double> arr = { datacopy.table[0].data[i],dQ };
		temp.append(arr);
	}
	return temp;
}

bool CylinderPressure::analyze(CylinderGeometry Geo)
{
	return true;
}

ArrayTable CylinderPressure::ploytropicIndex(CylinderGeometry CylGeo)
{
	ArrayTable result(2, 0);
	for (int i = 0; i < data.row - 1; i++)
	{
		std::vector<double> temp = { data.table[0].data[i],log(data.table[1].data[i + 1] / data.table[1].data[i]) / log(CylGeo.V(data.table[0].data[i]) / CylGeo.V(data.table[0].data[i + 1])) };
		result.append(temp);
	}
	return result;
}

ArrayTable CylinderPressure::slice2CombustionPeriod(double _intakeValeClose, double _exhaustValveOpen)
{
	int maxCA = data.findMaxDataIndex(1);
	_intakeValeClose = _intakeValeClose > data.table[0].data[maxCA] ? _intakeValeClose - 720 : _intakeValeClose;
	_exhaustValveOpen = _exhaustValveOpen < data.table[0].data[maxCA] ? _exhaustValveOpen + 720 : _exhaustValveOpen;
	if (data.table[0].data[0] > _intakeValeClose)
	{
		for (int i = data.row - 1; data.table[0].data[i] >= _intakeValeClose + 720; i--)
		{
			data.table[0].data[i] = data.table[0].data[i] - 720;
		}
	}
	else
	{
		for (int i = 0; data.table[0].data[i] < _intakeValeClose; i++)
		{
			data.table[0].data[i] = data.table[0].data[i] + 720;
		}
	}
	data.doQuickSort(0);
	for (int i = data.row - 1; data.table[0].data[i] > _exhaustValveOpen; i--)
	{
		data.pop();
	}
	return data;
}

CylinderPressure::CylinderPressure(ArrayTable data) :data(data)
{}

double CylinderPressure::MultipointSmooth(int j, smoothtype type, ArrayTable input)
{
	double p[41];
	if (type == smoothtype::_10pointS)
	{
		for (int i = 0; i <= 10; i++)
		{
			int k = j + i - 5;
			if (k < 0)
			{
				k = input.row + k;
			}
			if (k >= input.row)
			{
				k = k - input.row;
			}
			p[i] = input.table[1].data[k];
		}
		return 0.01406228988022 * (p[0] + p[10]) + 0.02999161403147 * (p[1] + p[9])
			+ 0.07199291837308 * (p[2] + p[8]) + 0.12455311765726 * (p[3] + p[7])
			+ 0.16744779125392 * (p[4] + p[6]) + 0.18390453760809 * p[5];
	}
	if (type == smoothtype::_20pointS)
	{
		for (int i = 0; i <= 20; i++)
		{
			int k = j + i - 10;
			if (k < 0)
			{
				k = input.row + k;
			}
			if (k >= input.row)
			{
				k = k - input.row;
			}
			p[i] = input.table[1].data[k];
		}
		return 0.00628875982703 * (p[0] + p[20]) + 0.00835789480574 * (p[1] + p[19])
			+ 0.01413216724080 * (p[2] + p[18]) + 0.02334549711773 * (p[3] + p[17])
			+ 0.03528665004430 * (p[4] + p[15]) + 0.04886897162590 * (p[5] + p[15])
			+ 0.06275158799024 * (p[6] + p[14]) + 0.07549718884042 * (p[7] + p[13])
			+ 0.08574616653274 * (p[8] + p[12]) + 0.09238429006941 * (p[9] + p[11])
			+ 0.09468165181140 * p[10];
	}
	if (type == smoothtype::_30pointS)
	{
		for (int i = 0; i <= 30; i++)
		{
			int k = j + i - 15;
			if (k < 0)
			{
				k = input.row + k;
			}
			if (k >= input.row)
			{
				k = k - input.row;
			}
			p[i] = input.table[1].data[k];
		}
		return 0.00337275314585 * (p[0] + p[30]) + 0.00403927236260 * (p[1] + p[29])
			+ 0.00568295044384 * (p[2] + p[28]) + 0.00837876896325 * (p[3] + p[27])
			+ 0.01213324475635 * (p[4] + p[26]) + 0.01687617855219 * (p[5] + p[25])
			+ 0.02246326348381 * (p[6] + p[24]) + 0.02868209834391 * (p[7] + p[23])
			+ 0.03526379680823 * (p[8] + p[22]) + 0.04189908850144 * (p[9] + p[21])
			+ 0.04825780466797 * (p[10] + p[20]) + 0.05401036038131 * (p[11] + p[19])
			+ 0.05884967029500 * (p[12] + p[18]) + 0.06251188108049 * (p[13] + p[17])
			+ 0.06479437631446 * (p[14] + p[15]) + 0.06556970379859 * p[15];
	}
	if (type == smoothtype::_40pointS)
	{
		for (int i = 0; i <= 40; i++)
		{
			int k = j + i - 20;
			if (k < 0)
			{
				k = input.row + k;
			}
			if (k >= input.row)
			{
				k = k - input.row;
			}
			p[i] = input.table[1].data[k];
		}
		return 0.00174740765868 * (p[0] + p[40]) + 0.00207106934517 * (p[1] + p[39])
			+ 0.00271625194969 * (p[2] + p[38]) + 0.00374632977076 * (p[3] + p[37])
			+ 0.00520958015585 * (p[4] + p[36]) + 0.00713555546439 * (p[5] + p[35])
			+ 0.00953230945972 * (p[6] + p[34]) + 0.01238465536383 * (p[7] + p[33])
			+ 0.01565358008326 * (p[8] + p[32]) + 0.01927687844556 * (p[9] + p[31])
			+ 0.02317100599510 * (p[10] + p[30]) + 0.02723408261422 * (p[11] + p[29])
			+ 0.03134991568708 * (p[12] + p[28]) + 0.03539285432105 * (p[13] + p[27])
			+ 0.03923323856728 * (p[14] + p[26]) + 0.04274317240413 * (p[15] + p[25])
			+ 0.04580232854741 * (p[15] + p[24]) + 0.04830348820519 * (p[17] + p[23])
			+ 0.05015753009378 * (p[18] + p[22]) + 0.05129760985250 * (p[19] + p[21])
			+ 0.05168231203067 * p[20];
	}
}


double BMEP(double n, double Pe, double TVH, int stroke)
{
	return Pe * (30 * stroke / n) / TVH;
}


double FMEP(double D, double cm, double pme)
{
	return pow(D, -0.1778) * (0.00855 * cm + 0.0789 * pme/1.e6 - 0.0214)*1.e6;//MPa
}

double MeEff(double D, double cm, double pme)
{
	double fmep = FMEP(D, cm, pme);
	return pme / (pme + fmep);
}

double PIntakeManifold(double Ts, double gf, double VC, double alpha, double VolEff,double L0)
{
	return justi::Rg(1.e8)* Ts* (L0 * gf * alpha) / VC / VolEff;
}

double BSFC(double eta_i, double eta_m, double Hu)
{
	return 3.6e6 / (Hu / 1.e3 * eta_i * eta_m);
}

ArrayTable BSFCExample()
{
	ArrayTable result(2, 0);
	result.setParaName(std::vector<std::string> {"Speed ratio,Speed/rate-speed", "Brake spacific fuel consumption"});
	result.setParaUnit(std::vector<std::string> {"/", "g/(kW·h)"});
	auto fun = [](double x)->double {return -827.92 * pow(x, 4) + 2163.5 * pow(x, 3) - 1835.2 * pow(x, 2) + 526.25 * x + 177.48; };
	for (double i = 0.2; i < 1.2; i+=0.01)
	{
		result.append(std::vector<double> {i, fun(i)});
	}
	return result;
}

ArrayTable IdealModel(CylinderGeometry Geo, double Rp, double alpha, double Ps, double Ts, double L0 , double Hu)
{
	ArrayTable result(3, 0); 
	result.setParaName(std::vector<std::string> {"Volumn", "Pressure", "Temperature"});
	result.setParaUnit(std::vector<std::string>{"m^3", "Pa", "K"});
	result.append(std::vector<double> {Geo.V(180), Ps, Ts});

	//等熵压缩阶段
	double Ttemp = Ts; double ptemp = Ps; double k = justi::k_Justi(Ttemp);
	for (double i = 181; i <= 360; i++)
	{
		k = justi::k_Justi(Ttemp,alpha);
		ptemp = ptemp * pow(Geo.V(i-1), k) / pow(Geo.V(i), k);
		Ttemp = Ttemp * pow(Geo.V(i-1), k - 1) / pow(Geo.V(i), k - 1);
		result.append(std::vector<double> {Geo.V(i), ptemp, Ttemp});
	}

	auto fun = [=](double Tz)->double {return (justi::cv_mean_Justi(Ttemp, Tz, alpha) * (1. + alpha * L0) * (Tz - Ttemp) - Rp * Hu)/1.e3; };
	/*ArrayTable te(2, 0);
	for (double i = Ttemp+1; i < 2000; i++)
	{
		te.append(std::vector<double> {i, fun(i)});
	}
	te.openWithProgram();*/
	double Tzinit = Ttemp + 100;
	while (abs(fun(Tzinit)) > 1.)
	{
		double h = 1.;
		double dfun = (fun(Tzinit + h) - fun(Tzinit - h)) / (2 * h);
		Tzinit = Tzinit - fun(Tzinit) / dfun;
	}
	double pinit = ptemp * Tzinit / Ttemp;
	result.append(std::vector<double> {Geo.V(360), pinit, Tzinit});

	double ma = Ps * Geo.V(180) / justi::Rg(1.e8) / Ts;

	auto fun2 = [=](double Tz)->double {return (justi::cp_mean_Justi(Tzinit, Tz, alpha) * (1. + 1/alpha / L0)*ma * (Tz - Tzinit)+pinit*Geo.V(360)*(1-Tz/Tzinit) - (1-Rp) * Hu* 1 / alpha / L0*ma) / 1.e3; };
	double Tzend = Tzinit + 100;
	while (abs(fun(Tzend)) > 1.)
	{
		double h = 1.;
		double dfun = (fun2(Tzend + h) - fun2(Tzend - h)) / (2 * h);
		Tzend = Tzend - fun(Tzend) / dfun;
		std::cout << "Tzend" << Tzend;
		
	}
	result.append(std::vector<double> {Geo.V(360)* Tzinit/Tzend, pinit, Tzend});
	

	return result;
}

heatRelease::Combustion::Combustion()
{
	data = ArrayTable(3, 0);
	data.setParaName(std::vector<std::string> {"cranke angle", "DX", "X"});
	data.setParaUnit(std::vector<std::string> {"CA", "/", "/"});
}
