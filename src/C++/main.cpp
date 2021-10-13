#include<pybind11/pybind11.h>
#include<pybind11/stl.h>

#include "Valve.h"

# include "Cylinder.h"

#include "Functions.h"

namespace py=pybind11;

using namespace py::literals;

PYBIND11_MODULE(ICELib,m){
    py::class_<CylinderGeometry>(m, "CylinderGeometry")
        .def(py::init<double, double, double, double>())
        .def("V", &CylinderGeometry::V)
        .def("DV",&CylinderGeometry::DV)
        .def("poltVolme",&CylinderGeometry::plotVolume)
        .def_readwrite("bore",&CylinderGeometry::D)
        .def_readwrite("stroke",&CylinderGeometry::S)
        .def_readwrite("compression ratio",&CylinderGeometry::Eps)
        .def_readwrite("data", &CylinderGeometry::data);

    py::class_<heatRelease::SingleWiebe>(m, "SingleWiebe")
        .def(py::init<double, double, double, double>())
        .def("DX",&heatRelease::SingleWiebe::DX)
        .def("X", &heatRelease::SingleWiebe::X)
        .def("getFractionCA", &heatRelease::SingleWiebe::getFractionCA)
        .def_readwrite("data",&heatRelease::SingleWiebe::data);

    py::class_<ArrayTable>(m, "ArrayTable")
        .def(py::init<int,int,std::string>())
        .def("show", py::overload_cast<>(&ArrayTable::show))
        .def("openWithProgram",&ArrayTable::openWithProgram)
        .def("plot",&ArrayTable::plot,"_coly"_a,"_colx"_a=0);

    ////ÈÝ»ýÀà
    py::class_<Volume>(m, "Volume")
        .def(py::init<double, double, double, double>())
        .def("showProp", &Volume::showProp)
        .def("DMDT", &Volume::DMDT)
        .def("DTDT", &Volume::DTDT, "DQDT"_a = 0)
        .def("DAFA", &Volume::DAFA)
        .def("M_exhaust", &Volume::M_exhasut)
        .def("recordThisStep", &Volume::recordThisStep)
        .def("updateStatus", &Volume::updateStatus)
        .def_readwrite("T", &Volume::T)
        .def_readwrite("m", &Volume::M)
        .def_readwrite("AFA", &Volume::AFA)
        .def_readwrite("p", &Volume::p)
        .def_readwrite("data",&Volume::data)
        .def("__repr__", [](const Volume& a) {
        return "A Volume whose volume is " + function::double2string(a.V)+"\n"
            +"temperature "+function::double2string(a.T)+" K\n"+
            "pressure "+function::double2string(a.p/1.e5)+" bar\n"; });

    /*py::class_<valve>(m,"valve")
        .def(py)*/

    py::class_<ValveSimple>(m, "ValveSimple")
        .def(py::init<double>())
        .def("massFlow",&ValveSimple::massFlow)
        .def("connect",&ValveSimple::connect);

    /*py::class_<volComponent>(m, "volComponent")
        .def(py::init<>());*/


    m.def("FMEP", &FMEP, "calculate the FMEP");

    m.def("FlowUnitArea", &flowUnitArea,"nozzle equation");

    m.def("ondas_out", &ondas_out, "ondas_out lamda_in,AA,phi,pp,k,");

    m.def("ondas_valve", &ondas_valve, "ondas_valve");

}

int main()
{
    ValveSimple* V =new ValveSimple();
    Volume* V1 = new Volume();
    Volume* V2 = new Volume();
    V->connect(V1, V2);
}

