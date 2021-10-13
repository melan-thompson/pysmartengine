#include "Valve.h"


ValveSimple::ValveSimple(double _flowArea) :flowArea(_flowArea)
{

}


double flowUnitArea(double P1, double T1, double R1, double K1, double P2, double T2, double R2, double K2)
{
	if (P1 < 0 || P2 < 0)
	{
		std::cerr << "Pressure can not less than 0!!!\n"; system("pause");
	}
	if (T1 < 0 || T2 < 0)
	{
		std::cerr << "Temperature can not less than 0!!!\n"; system("pause");
	}
	if (R1 < 0 || R2 < 0)
	{
		std::cerr << "Gas constant can not less than 0!!!\n"; system("pause");
	}
	if (K1 < 0 || K2 < 0)
	{
		std::cerr << "Inapropriate spacific heat ratio!!!\n"; system("pause");
	}

	double TMP2 = P2 / P1;

	//����Ϊ��
	if (TMP2 < 1.0)
	{
		double TMP1 = sqrt(2 * K1 / R1 / T1);

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

	//����Ϊ��
	else
	{
		double TMP1 = sqrt(2 * K2 / R2 / T2);
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

ArrayTable ValveSimple::airFlowExample(double P1, double T1, double P2, double T2)
{
	DieselMixture air(0.1, 0, 0.001);
	ArrayTable temp(2, 0);
	std::vector<std::string> name = { "Pressure ratio","Mass flow rate" };
	std::vector<std::string> unit = { "/","kg/s" };
	temp.setParaName(name); temp.setParaUnit(unit);
	temp.show();
	/*for (double i = 0; i < P1; i += 0.001 * P1)
	{
		std::vector<double> arr = { i / P1,flowArea * flowUnitArea(P1,T1,air.Rg(0),air.k(0,T1),i) };
		temp.append(arr);
	}*/
	//temp.openWithProgram();

	if (P2 != NULL)
	{
		for (double i = 0.3 * P2; i < P2; i += 0.001 * P2)
		{
			std::vector<double> arr = { P2 / i,-flowArea * flowUnitArea(i,T1,air.Rg(0),air.k(0,T1),P2,T2,air.Rg(0),air.k(0,T2)) };
			temp.append(arr);
		}

		temp.doQuickSort(0);
	}

	return temp;
}

double ValveSimple::massFlowRate(double P1, double T1, double R1, double K1, double P2, double T2, double R2, double K2)
{
	return flowArea * flowUnitArea(P1, T1, R1, K1, P2, T2, R2, K2);
}

double ValveSimple::massFlow()
{
	if (comLast==nullptr||comNext==nullptr)
	{
		std::cerr << "You must connect the valve to some components\n";
	}

	/*std::cout << comLast->get_p()<<std::endl;
	std::cout << comNext->get_p() << std::endl;*/

	DMDT = flowArea * flowUnitArea(comLast->get_p(), comLast->get_T(), comLast->get_R(), comLast->get_k(), comNext->get_p(), comNext->get_T(), comNext->get_R(), comNext->get_k());
	return DMDT;
}


bool ValveSimple::connect(volComponent* connectFrom, volComponent* connectTo)
{
	volComponent* connectFrom1 = connectFrom;

	volComponent* connectTo1 = (volComponent*) connectTo;

	//上游的下游节点为本节点,连接的上游节点为

	//flowLink* temp = (flowLink*)this;

	connectFrom->flowNext.push_back(this);
	comLast = connectFrom1;

	//下游的上游节点为本节点
	connectTo->flowLast.push_back(this);
	comNext = connectTo1;

	std::cout << "Connect successfully\n";
	return true;
}



Valve::Valve(ArrayTable Valve_LiftIn, double valveDiameter, double Open_Timing_Angle, double Close_Timing_Angle, double SigmaIn, double rock, ArrayTable _flow_coeff)
	: valveDiameter(valveDiameter), Open_Timing_Angle(Open_Timing_Angle), Close_Timing_Angle(Close_Timing_Angle),
	rock(rock), Valve_Lift(Valve_LiftIn), flowCoeff(_flow_coeff)
{
	//��ʼ����������
	double open = Valve_LiftIn.table[0].data[0], close = Valve_LiftIn.table[0].data[Valve_LiftIn.row - 1];
	double temp = (Close_Timing_Angle - Open_Timing_Angle) / (close - open);
	for (int i = 0; i < Valve_LiftIn.row; i++)
	{
		Valve_Lift.table[0].data[i] = Open_Timing_Angle + temp * (Valve_LiftIn.table[0].data[i] - open);
		Valve_Lift.table[1].data[i] = Valve_LiftIn.table[1].data[i] * rock;
	}

	sigma = SigmaIn * PI / 180;

	record = ArrayTable(3, 0);
	std::vector < std::string>  paraName = { "����ת��","��������","��������" };
	std::vector < std::string> paraUnit = { "CA","m","kg/s" };
	record.setParaName(paraName);
	record.setParaUnit(paraUnit);

	flowCoeff.setParaName(std::vector<std::string> {"Crank angle", "Flow coefficiency"});
	flowCoeff.setParaUnit(std::vector<std::string> {"CA", "/"});

	
	double maxlift = Valve_Lift.findMaxData(1);
	flowCoeff.Clear();
	for (double i = 0; i <= maxlift; i+=0.001)
	{
		flowCoeff.append(std::vector<double> {i, flowcoefficient(i)});
	}
	

};

//��������ϵ���ľ��鹫ʽ
double Valve::flowcoefficient(double _valveLift)
{
	return 0.98 - 3.3 * pow((_valveLift / valveDiameter), 2);
}
