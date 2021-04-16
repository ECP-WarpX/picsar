#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "picsar_qed/physics/chi_functions.hpp"

#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_core.hpp"
#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"

#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_core.hpp"
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"

#include "picsar_qed/physics/schwinger/schwinger_pair_engine_core.hpp"

#include <vector>

namespace py = pybind11;

namespace pxr_phys = picsar::multi_physics::phys;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;
namespace pxr_sc = picsar::multi_physics::phys::schwinger;

PYBIND11_MODULE(pxr_qed, m) {
    m.doc() = "pybind11 pxr_qed plugin";

#ifdef PXRQEDPY_DOUBLE_PRECISION
    using PXRQEDPY_REAL = double;
    m.attr("PRECISION") = py::str("double");
#else
    using PXRQEDPY_REAL = float;
    m.attr("PRECISION") = py::str("float");
#endif

#if PXRQEDPY_UNITS == SI
    constexpr auto PXRQEDPY_UU = pxr_phys::unit_system::SI;
    m.attr("UNITS") = py::str("SI");
#elif PXRQEDPY_UNITS == NORM_OMEGA
    constexpr auto PXRQEDPY_UU = pxr_phys::unit_system::norm_omega;
    m.attr("UNITS") = py::str("NORM_OMEGA");
#elif PXRQEDPY_UNITS == NORM_LAMBDA
    constexpr auto PXRQEDPY_UU = pxr_phys::unit_system::norm_lambda;
    m.attr("UNITS") = py::str("NORM_LAMBDA");
#elif PXRQEDPY_UNITS == HEAVISIDE_LORENTZ
    constexpr auto PXRQEDPY_UU = pxr_phys::unit_system::heaviside_lorentz;
    m.attr("UNITS") = py::str("HEAVISIDE_LORENTZ");
#else
    #error PXRQEDPY_UNITS is incorrectly defined!
#endif

    using VECTOR = std::vector<PXRQEDPY_REAL>;
    using RAW_DATA = std::vector<char>;

// ******************************* Chi parameters ***********************************************
    m.def("chi_photon",
        py::vectorize(       
            py::overload_cast<
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL>(
                &pxr_phys::chi_photon<PXRQEDPY_REAL, PXRQEDPY_UU>)),
        "Returns the chi parameter for a photon");

    m.def("chi_ele_pos",
        py::vectorize(       
            py::overload_cast<
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL>(
                &pxr_phys::chi_ele_pos<PXRQEDPY_REAL, PXRQEDPY_UU>)),
        "Returns the chi parameter for an electron or a positron");
// ______________________________________________________________________________________________



// ******************************* Breit-Wheeler pair production ********************************
    auto bw = m.def_submodule( "bw" );

    using bw_dndt_lookup_table_params =
        pxr_bw::dndt_lookup_table_params<PXRQEDPY_REAL>;

    py::class_<bw_dndt_lookup_table_params>(bw,
        "dndt_lookup_table_params",
        "Parameters to generate a dN/dt lookup table")
        .def(py::init<>())
        .def(py::init<PXRQEDPY_REAL,PXRQEDPY_REAL,int>())
        .def("__eq__", &bw_dndt_lookup_table_params::operator==)
        .def_readwrite("chi_phot_min", &bw_dndt_lookup_table_params::chi_phot_min)
        .def_readwrite("chi_phot_max", &bw_dndt_lookup_table_params::chi_phot_max)
        .def_readwrite("chi_phot_how_many", &bw_dndt_lookup_table_params::chi_phot_how_many);

    using bw_dndt_lookup_table =
        pxr_bw::dndt_lookup_table<PXRQEDPY_REAL, VECTOR>;

    py::class_<bw_dndt_lookup_table>(bw,
        "dndt_lookup_table",
        "dN/dt lookup table")
        .def(py::init<bw_dndt_lookup_table_params>())
        .def("__eq__", &bw_dndt_lookup_table::operator==)
        .def("generate", &bw_dndt_lookup_table::generate)
        .def("serialize", &bw_dndt_lookup_table::serialize);    

    using bw_pair_prod_lookup_table_params =
        pxr_bw::pair_prod_lookup_table_params<PXRQEDPY_REAL>;

    py::class_<bw_pair_prod_lookup_table_params>(bw, "pair_prod_lookup_table_params")
        .def(py::init<>())
        .def(py::init<PXRQEDPY_REAL,PXRQEDPY_REAL,int,int>())
        .def("__eq__", &bw_pair_prod_lookup_table_params::operator==)
        .def_readwrite("chi_phot_min", &bw_pair_prod_lookup_table_params::chi_phot_min)
        .def_readwrite("chi_phot_max", &bw_pair_prod_lookup_table_params::chi_phot_max)
        .def_readwrite("chi_phot_how_many", &bw_pair_prod_lookup_table_params::chi_phot_how_many)
        .def_readwrite("frac_how_many", &bw_pair_prod_lookup_table_params::frac_how_many);

    bw.def("get_optical_depth", 
        py::vectorize(pxr_bw::get_optical_depth<PXRQEDPY_REAL>), 
        "Computes the optical depth of a new photon");

    bw.def("get_dn_dt", 
            [](std::vector<PXRQEDPY_REAL> energy_phot, std::vector<PXRQEDPY_REAL> chi_phot,
                const bw_dndt_lookup_table& ref_table, PXRQEDPY_REAL ref_quantity){
                    assert(energy_phot.size() == chi_phot.size());
                    assert(ref_table.is_init());

                    auto res = std::vector<PXRQEDPY_REAL>(energy_phot.size());

                    std::transform(
                        energy_phot.begin(),
                        energy_phot.end(), 
                        chi_phot.begin(), 
                        res.begin(), 
                        [&](PXRQEDPY_REAL ee, PXRQEDPY_REAL chi){
                            return pxr_bw::get_dN_dt<
                                PXRQEDPY_REAL, bw_dndt_lookup_table, PXRQEDPY_UU>(
                                ee, chi,
                                ref_table, ref_quantity,nullptr);
                        });

                    return res;                   
                ;}, 
        "Computes dN/dt for Breit-Wheeler pair production");
// ______________________________________________________________________________________________



// ******************************* Quantum syncrhotron emission *********************************
    auto qs = m.def_submodule( "qs" );

    qs.def("get_optical_depth", 
        py::vectorize(pxr_qs::get_optical_depth<PXRQEDPY_REAL>), 
        "Computes the optical depth of a new electron or positron");

// ______________________________________________________________________________________________



// ******************************* Schwinger pair production *************************************
    auto sc = m.def_submodule( "sc" );

    sc.def("pair_production_rate",
        py::vectorize(       
            py::overload_cast<
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL>(
                    &pxr_sc::pair_production_rate<PXRQEDPY_REAL, PXRQEDPY_UU>)),
        "Computes the Schwinger pair production rate using the Nikishov formula");

    sc.def("expected_pair_number",
        py::vectorize(       
            py::overload_cast<
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL,
                PXRQEDPY_REAL,PXRQEDPY_REAL,PXRQEDPY_REAL>(
                    &pxr_sc::expected_pair_number<PXRQEDPY_REAL, PXRQEDPY_UU>)),
        "Computes the Schwinger pair production rate using the Nikishov formula");
// ______________________________________________________________________________________________

}
