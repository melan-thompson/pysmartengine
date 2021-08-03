#include <pybind11/pybind11.h>

namespace py = pybind11;


double u_Justi(double T, double AFAK, double Tref)
{
	double temp1 = 489.6 + 46.4 / pow(AFAK, 0.93);
	double temp2 = 7.768 + 3.36 / pow(AFAK, 0.8);
	double temp3 = 0.0975 + 0.0485 / pow(AFAK, 0.75);
	double temp4 = 1356.8 + temp1 * (T - Tref) * 1.e-2 + temp2 * pow(T - Tref, 2) * 1.e-4 - temp3 * pow(T - Tref, 3) * 1.e-6;
	return 0.1445 * temp4 * 1.e3;//J/kg
}

double Rg(double AFAK)
{
	return 9.81 * (29.2647 - 0.0402 / AFAK);//J/(kg*K)
}


PYBIND11_MODULE(pybind2, m) {
	m.doc() = "calculating fluid properties";
	m.def("u_Justi", &u_Justi, "u_Justi");
	//m.def("Rg", &Rg, "Rg");
}