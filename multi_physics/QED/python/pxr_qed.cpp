#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_core.hpp"
#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"

#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_core.hpp"
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"

#include "picsar_qed/physics/schwinger/schwinger_pair_engine_core.hpp"

#define PXRQEDPY_REAL double
#define PXRQEDPY_UNITS pxr_phys::unit_system::SI

namespace py = pybind11;

namespace pxr_phys = picsar::multi_physics::phys;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;
namespace pxr_sc = picsar::multi_physics::phys::schwinger;

PYBIND11_MODULE(pxr_qed, m) {
    m.doc() = "pybind11 pxr_qed plugin";

    using REAL = PXRQEDPY_REAL;

// ******************************* Breit-Wheeler pair production ********************************
    auto bw = m.def_submodule( "bw" );

    using bw_dndt_lookup_table_params =
        pxr_bw::dndt_lookup_table_params<REAL>;

        py::class_<bw_dndt_lookup_table_params>(bw, "dndt_lookup_table_params")
            .def(py::init<>())
            .def(py::init<REAL,REAL,int>())
            .def("__eq__", &bw_dndt_lookup_table_params::operator==)
            .def_readwrite("chi_phot_min", &bw_dndt_lookup_table_params::chi_phot_min)
            .def_readwrite("chi_phot_max", &bw_dndt_lookup_table_params::chi_phot_max)
            .def_readwrite("chi_phot_how_many", &bw_dndt_lookup_table_params::chi_phot_how_many);
       
    bw.def("get_optical_depth", 
        py::vectorize(pxr_bw::get_optical_depth<REAL>), 
        "Computes the optical depth of a new photon");

// ______________________________________________________________________________________________


// ******************************* Quantum syncrhotron emission *********************************
    auto qs = m.def_submodule( "qs" );

    qs.def("get_optical_depth", 
        py::vectorize(pxr_qs::get_optical_depth<REAL>), 
        "Computes the optical depth of a new electron or positron");

// ______________________________________________________________________________________________



// ******************************* Schwinger pair production *************************************
    auto sc = m.def_submodule( "sc" );


    sc.def("pair_production_rate",
        py::vectorize(       
            py::overload_cast<REAL,REAL,REAL,REAL,REAL,REAL, REAL>(
                &pxr_sc::pair_production_rate<REAL, PXRQEDPY_UNITS>)),
        "Computes the Schwinger pair production rate using the Nikishov formula");

    sc.def("expected_pair_number",
        py::vectorize(       
            py::overload_cast<REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL>(
                &pxr_sc::expected_pair_number<REAL, PXRQEDPY_UNITS>)),
        "Computes the Schwinger pair production rate using the Nikishov formula");
// ______________________________________________________________________________________________

}
