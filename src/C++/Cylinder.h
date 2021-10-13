#ifndef CYLINDER_H
#define CYLINDER_H
#include<cmath>
#include "ArrayTable.h"
#include "ParameterTable.h"
#include "GasProperty.h"
#include "FlowComponent.h"

double ondas_out(double lambda_in, double AA, double phi, double rb, double k, double tor);

/// <summary>
/// 流入管道计算Aa,lambdain,lambdaout
/// </summary>
/// <param name="clin"></param>
/// <param name="clout"></param>
/// <param name="AAn"></param>
/// <param name="phi"></param>
/// <param name="pp">p/pref</param>
/// <param name="A0">外界声速</param>
/// <param name="k"></param>
/// <param name="INFLOW_CONST_P"></param>
/// <param name="tolLoop"></param>
/// <returns></returns>
std::vector<double> ondas_valve(double clin, double clout, double AAn, double phi, double pp, double A0, double k, bool INFLOW_CONST_P, double tolLoop);


namespace heatRelease 
{
	class Combustion
	{
	public:
		ArrayTable data;
		double SOC, EOC;
		Combustion();
		virtual double X(double CA) = 0;
		virtual double DX(double CA) = 0;
		virtual double getFractionCA(double fraction) = 0;
	};

	/// <summary>
	/// <param name="CA">曲轴转角</param>
	/// <param name="SOC">Start of combustion</param>
	/// <param name="TCD">Total combustion duration</param>
	/// <param name="m">m,adjustable parameters</param>
	/// <param name="a">a is usually 6.908 by default</param>
	/// </summary>
	class SingleWiebe:public Combustion
	{
	public:
		SingleWiebe() {};

		SingleWiebe(double SOC, double TCD, double m, double a = 6.908);

		double SOC, TCD, m, a;

		double X(double CA);

		double DX(double CA);

		double getFractionCA(double fraction=0.5f);

	};

	class DoubleWiebe:public Combustion
	{
	public:
		DoubleWiebe(double PF, double PSOC, double PTCD, double DTCD, double PM = 2, double DM = 0.8, double DSOCin = INFINITY);

		double PF;//Premixed combustion fraction
		double DF;//Disfusion combustion fraction
		double PSOC;//Premixed combustion start of combustion
		double PTCD;//Premixed combustion total combustion duration
		double PM;//Premixed combusiton exponent
		double DSOC;//Disfusion combustion start of combustion 
		double DTCD;//Disfusion combustion total combustion duration
		double DM;//Disfusion combustion exponent

		SingleWiebe PremixWiebe;
		SingleWiebe DisfusionWiebe;

		double X(double CA);

		double DX(double CA);

		double getFractionCA(double fraction = 0.5f);

	};

	class MultiWiebe :public Combustion
	{
	public:
		MultiWiebe(std::vector<double> SOC,std::vector<double> dur,std::vector<double> m);

		double X(double CA);

		double DX(double CA);

		double getFractionCA(double fraction = 0.5f);


	};


	class HeatReleaseData//气缸的燃烧放热模型
	{
	public:
		ArrayTable heat_release_data;
		HeatReleaseData() {};
		HeatReleaseData(ArrayTable data);
		~HeatReleaseData() {};

		double SOC();//Start of combustion
		

		double EOC();//End of combustion
			

		double DX(double Fi);
		

		double X(double Fi);
		
		//void fitWithSingleWiebe()
		//{
		//	/*static ArrayTable error(2, 0),data=heat_release_data;
		//	static double SOCIndex = findSOCIndex(), EOCIndex = findEOCIndex();
		//	static double _SOC = SOC(), _EOC = EOC();
		//	auto fun = [](std::vector<double> x)->double {
		//		SingleWiebe temp(_SOC, _EOC - _SOC, x[0]);
		//		for (int i = SOCIndex; i < EOCIndex; i++)
		//		{
		//			std::vector<double> arr = { data.table[0].data[i],
		//				data.table[1].data[i] - temp.DX(data.table[0].data[i]) };
		//			error.append(arr);
		//		};
		//		return error.quasierror(1);
		//	}; 

		//	
		//	std::vector<double> lb = { 0.1 };
		//	std::vector<double> up = { 5 };
		//	
		//	ParticleGroup train(fun, 100, lb, up);
		//	train.go();*/
		//

		//}

		/*double fun(std::vector<double> x)
		{
			SingleWiebe temp(SOC(), EOC() - SOC(), x[0]);
			ArrayTable error(2, 0);
			for (int i = findSOCIndex(), int end = findEOCIndex(); i < end; i++)
			{
				std::vector<double> arr = { heat_release_data.table[0].data[i],
					heat_release_data.table[1].data[i] - temp.DX(heat_release_data.table[0].data[i]) };
				error.append(arr);
			}
			return error.quasierror(1);
		}*/

		/*std::string fitWithDoubleWiebe()
		{

		}*/

		ParameterTable analyze();
	private:
		int findSOCIndex();//Start of combustion

		int findEOCIndex();//End of combustion

	};

}

namespace heatTransfer
{
	class CylinderHeatTransfer//气缸的传热模型
	{
	public:

		CylinderHeatTransfer() {};
		CylinderHeatTransfer(double Head_Temperature, double Piston_Temperature, double Cylinder_Temperature, double CorrectCoeff = 1);

		//Pe为平均有效压力
		CylinderHeatTransfer(double Pe, double CorrectCoeff = 1);

		double HeatTrasferRate(double S, double D, double V, double N, double P, double T);

	private:
		double HT;// Head_Temperature;
		double PT;// Piston_Temperature;
		double CT;// Cylinder_Temperature;
		double CQW;//Heat transfer coefficient

	};

	///<summary>Sitkei传热公式,适用于小型柴油机计算</summary>
	double Sitkei(double p, double T, double Cm, double D, double S, double V);

	///<summary>Woschni公式</summary>
	double simpleWoschini(double D, double vm, double pz, double Tz);

}

class CylinderGeometry//气缸的几何参数模型，已验证其正确性
{
public:
	double D;//Bore;
	double S;//Stroke;
	double Eps;//Compression_Ratio;
	double Lam;// Lambda;
	ArrayTable data;

	CylinderGeometry() {};
	CylinderGeometry(double Bore, double Stroke, double CR, double length);
	~CylinderGeometry() {};
	double VH;

	double V(double Fi) const;

	double DV(double Fi) const;

	double clearanceVolume() const;

	double TDCclearanceHeight() const;

	double heatExchangeArea(double crank_angle)const;

	void plotVolume();
		
private:
	double STRC2, STRC3, STRC4, STRC5;
};

class Injection//气缸喷油模型
{
public:
	double MF;
	double L0;
	double QH;
	Injection() {};
	Injection(double MF, double L0, double QH);
};

class Volume:public volComponent//气缸的热物性
{
public:
	double V, T, M, AFA;//四个参数确定容积状态
	double p,rho;//pressure,密度
	double Rg, k, cp, cv;//gas property, all by calculating
	double h, u;//焓和热力学能
	ArrayTable data;
	Volume() {};

	//由容积压力，温度，体积，和过量空气系数来初始化容积状态量
	Volume(double Pressure, double Temperature, double Volume, double AFA);

	//更新状态参数
	bool updateStatus();

	void recordThisStep(double t);

	//bool update();

	//展示容积内的状态表
	void showProp();

	//计算质量的变化
	double DMDT();

	//计算温度的变化,首先要算质量的变化才能算温度的变化，默认传热量为0
	double DTDT(double DQDT=0);

	double DAFA();

	double M_exhasut();

	double M_air();

	

private:

	//时间步长
	double delta_time;

	double dm;//记录质量的变化

	double _Z(double p, double T, double AFAK);//已验证其正确性
	

	double _Rg(double AFAK);

	double cv_Justi(double T, double AFAK, double Tref = 273.15);

	double _M(double AFAK);
	


	double _k(double T, double AFAK);
	

	double cp_Justi(double T, double AFAK);
	

	double u_Justi(double T, double AFAK, double Tref = 273.15);

	double get_h() const { return h; };

	double get_AFA() const { return AFA; };

	double get_p()const { return p; };

	double get_T()const { return T; };

	double get_R()const { return Rg; };

	double get_k()const { return k; };

	const double L0=14.3;
	
};


/////////初始化就确定了的参数直接更新
class Cylinder :public Volume
{
public:
	double Fi;//曲轴转角local
	double N;//发动机转速
	double Xk;
	double DX;//算DM要用到，算DT要用到，算DAFAK要用到，避免重复计算，只计算一次
	ArrayTable record;
	Cylinder(double Fi, double pin, double Tin, double AFAin, double N, CylinderGeometry Geo, heatRelease::HeatReleaseData HRR, heatTransfer::CylinderHeatTransfer CylHT, Injection Inj);//初始化缸内参数，曲轴转角，压力，温度

	bool RecordThisStep();
	

	bool ClearRecord();
	
	std::vector<double> D(double DMI, double DMO, double DHI, double DHO, double AFAI, double AFAO);

	double M_exhaust();
	
	double M_air();
	

private:
	Injection Inj;
	CylinderGeometry Geo;
	heatRelease::HeatReleaseData HRR;
	heatTransfer::CylinderHeatTransfer CylHT;
	double _DMI_, _DMO_, _DHI_, _DHO_, _DQ_, _DQW_, _DW_;
	double _DM(double DMI, double DMO);

	double _DT(double DHI, double DHO, double DQW, double DQ, double DW, double DM);

	double _DMSC(double DMI, double DMO, double AFAOUT);

	double _DAFA(double DMI, double AFAI, double DMO, double AFAO, double DX);

	double _DQW();

};

class CylinderPressure
{
public:
	enum class smoothtype { _10pointS, _20pointS, _30pointS, _40pointS,_5point3times };
	ArrayTable data;
	ArrayTable dataAfterSmooth;
	CylinderPressure() {};
	CylinderPressure(ArrayTable data);
	bool readInFile(std::string filename);

	ArrayTable doMultipointSmooth(smoothtype type);

	double findStartOfCombustion();

	void plot();

	ArrayTable LogPLogV(CylinderGeometry CylGeo);

	ArrayTable PVDiagram(CylinderGeometry CylGeo);

	ArrayTable netHeatReleaseRate(CylinderGeometry CylGeo,DieselMixture mix);

	bool analyze(CylinderGeometry Geo);

	ArrayTable ploytropicIndex(CylinderGeometry CylGeo);

	ArrayTable slice2CombustionPeriod(double _intakeValeClose, double _exhaustValveOpen);

private:
		
	double MultipointSmooth(int j, smoothtype type, ArrayTable input);
};

/// <summary>
/// 计算平均有效压力
/// </summary>
/// <param name="n">Engine speed,发动机转速(r/min)</param>
/// <param name="Pe">Brake power,发动机功率(W)</param>
/// <param name="TVH">Total displacement volumn,发动机总排量(m^3)=PI/4*D*D*S</param>
/// <param name="stroke">Stroke,发动机冲程，默认为4</param>
/// <returns>平均有效压力(Pa)</returns>
double BMEP(double n,double Pe, double TVH,int stroke=4);

/// <summary>
/// 计算平均机械损失压力，适用于四冲程
/// </summary>
/// <param name="D">Bore，气缸直径(m)</param>
/// <param name="cm">Mean piston moving velocity，活塞平均运动速度(m/s)</param>
/// <param name="pme">Brake mean effective pressure，平均有效压力(Pa)</param>
/// <returns>FMEP，平均机械损失压力(Pa)</returns>
double FMEP(double D, double cm, double pme);

/// <summary>
/// 机械效率
/// </summary>
/// <param name="D">Bore，气缸直径(m)</param>
/// <param name="cm">Mean piston moving velocity，活塞平均运动速度(m/s)</param>
/// <param name="pme">Brake mean effective pressure，平均有效压力(Pa)</param>
/// <returns>Mechanical efficiency,机械效率=BMEP/(BMEP+FMEP)</returns>
double MeEff(double D, double cm, double pme);

/// <summary>
/// 估算进气压力
/// </summary>
/// <param name="Ts">Intake manifold temperature,进气管温度(K)</param>
/// <param name="gf">Inject fuel mass,单缸循环喷油量(kg)</param>
/// <param name="VC">Clearance volume,下止点容积(m^3)</param>
/// <param name="alpha">Excess air fuel ratio,过量空气系数</param>
/// <param name="VolEff">Volumetric efficiency,充量系数</param>
/// <returns>Intake manifold pressure,进气管压力(Pa)</returns>
double PIntakeManifold(double Ts, double gf, double VC, double alpha, double VolEff,double L0=14.3);


/// <summary>
/// 计算燃油消耗率
/// </summary>
/// <param name="eta_i">Indicated thermal efficiency,指示热效率</param>
/// <param name="eta_m">Mechanical efficiency,机械效率</param>
/// <param name="Hu">Low heat value,燃料低热值(J/kg)</param>
/// <returns>Brake spacific fuel consumption,有效燃油消耗率(g/(kW*h))</returns>
double BSFC(double eta_i, double eta_m, double Hu = 42496e3);

/// <summary>
/// 发动机油耗与发动机转速之间的关系
/// </summary>
/// <returns>1.转速比,发动机转速/标定转速，2.油耗(g/(kW*h))</returns>
ArrayTable BSFCExample();

ArrayTable IdealModel(CylinderGeometry Geo, double Rp, double alpha, double Ps, double Ts = 300, double L0 = 14.3, double Hu = 42496e3);
//{
//	ArrayTable result(3, 0);
//	result.setParaName(std::vector<std::string> {"Crank angle", "Pressure", "Temperature"});
//	result.setParaUnit(std::vector<std::string>{"CA", "Pa", "K"});
//	result.append(std::vector<double> {Geo.V(0), Ps, Ts});
//
//	//等熵压缩阶段
//	double Ttemp = Ts; double ptemp = Ps; double k = justi::k_Justi(Ttemp);
//	for (double i = 1; i <= 180; i++)
//	{
//		double ptemp = ptemp * pow(Geo.V(i), k) / pow(Geo.V(i - 1), k);
//		double Ttemp = Ttemp * pow(Geo.V(i), k - 1) / pow(Geo.V(i - 1), k - 1);
//		result.append(std::vector<double> {Geo.V(i), ptemp, Ttemp});
//	}
//
//	auto fun = [=](double Tz)->double{return justi::cv_mean_Justi(Ttemp, Tz,alpha)*(1. + alpha * L0)*(Tz - Ttemp) - Rp*Hu; };
//	double Tzinit = Ttemp + 30;
//	while (abs(Tzinit-Ttemp)>1.e-3)
//	{
//		double h = 1.e-2;
//		double dfun = (fun(Tzinit + h) - fun(Tzinit - h)) / (2 * h);
//		Tzinit = Tzinit - fun(Tzinit)/dfun;
//	}
//	result.append(std::vector<double> {Geo.V(180), ptemp* Tzinit / Ttemp, Tzinit});
//}
#endif