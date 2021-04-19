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
#include <string>
#include <sstream>
#include <tuple>
#include <iostream>
#include <fstream>

namespace py = pybind11;

namespace pxr_phys = picsar::multi_physics::phys;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;
namespace pxr_sc = picsar::multi_physics::phys::schwinger;

#define PXRQEDPY_DOUBLE_PRECISION 1

#ifdef PXRQEDPY_DOUBLE_PRECISION
    using REAL = double;
    const auto PXRQEDPY_PRECISION_STRING = std::string{"double"};
#else
    using REAL = float;
    const auto PXRQEDPY_PRECISION_STRING = std::string{"single"};
#endif

#ifdef PXRMP_HAS_OPENMP
    template <typename Func>
    void PXRQEDPY_FOR(int N, const Func& func){
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) func(i);
    }
    const auto PXRQEDPY_OPENMP_FLAG = true;
#else
    template <typename Func>
    void PXRQEDPY_FOR(int N, const Func& func){
        for (int i = 0; i < N; ++i) func(i);
    }
    const auto PXRQEDPY_OPENMP_FLAG = true;
#endif

#if PXRQEDPY_UNITS == SI
    const auto UU = pxr_phys::unit_system::SI;
    const auto PXRQEDPY_USTRING = std::string{"SI"};
#elif PXRQEDPY_UNITS == NORM_OMEGA
    const auto UU = pxr_phys::unit_system::norm_omega;
    const auto PXRQEDPY_USTRING = std::string{"NORM_OMEGA"};
#elif PXRQEDPY_UNITS == NORM_LAMBDA
    const auto UU = pxr_phys::unit_system::norm_lambda;
    const auto PXRQEDPY_USTRING = std::string{"NORM_LAMBDA"};
#elif PXRQEDPY_UNITS == HEAVISIDE_LORENTZ
    const auto UU = pxr_phys::unit_system::heaviside_lorentz;
    const auto PXRQEDPY_USTRING = std::string{"HEAVISIDE_LORENTZ"};
#else
    #error PXRQEDPY_UNITS is incorrectly defined!
#endif

template <typename Real>
inline std::string float_to_string(Real num)
{
    std::stringstream ss;
    ss << num;
    return ss.str();
}

inline std::string bool_to_string(bool val)
{
    return (val)?"True":"False";
}

using pyArr = py::array_t<REAL>;
using pyBufInfo = py::buffer_info;
using stdVec = std::vector<REAL>;
using rawVec = std::vector<char>;

void throw_error(const std::string& err_msg)
{
    throw std::runtime_error(" [ Error! ] " + err_msg);
}

template<typename ...Args>
auto aux_check_and_get_pointers(const long int len, const pyArr& last)
{
    const auto last_buf = last.request();
    if (last_buf.ndim != 1 || last_buf.shape[0] != len)
        throw_error("All arrays must be one-dimensional with equal size");

    const auto cptr = static_cast<REAL*>(last_buf.ptr);
    return std::make_tuple(cptr);
}

template<typename ...Args>
auto aux_check_and_get_pointers(const long int len, const pyArr& arg, const Args& ...args)
{
    const auto arg_buf = arg.request();
    if (arg_buf.ndim != 1 || arg_buf.shape[0] != len)
        throw_error("All arrays must be one-dimensional with equal size");

    const auto cptr = static_cast<REAL*>(arg_buf.ptr);
    return std::tuple_cat(std::make_tuple(cptr), 
        aux_check_and_get_pointers(len, args...));
}

template<typename ...Args>
auto check_and_get_pointers(const pyArr& first, const Args& ...args)
{
    const auto first_buf = first.request();
    if (first_buf.ndim != 1)
        throw_error("All arrays must be one-dimensional with equal size");
    
    const auto len = first_buf.shape[0];

    const auto cptr = static_cast<REAL*>(first_buf.ptr);

    return std::tuple_cat(std::make_tuple(len, cptr),
            aux_check_and_get_pointers(len, args...));
}

auto check_and_get_pointers(const pyArr& arr)
{
    const auto arr_buf = arr.request();
    if (arr_buf.ndim != 1)
        throw_error("Array must be one-dimensional");
    
    const auto len = arr_buf.shape[0];

    const auto cptr = static_cast<REAL*>(arr_buf.ptr);

    return std::make_tuple(len, cptr);
}

auto check_and_get_pointer_nonconst(pyArr& arr, const long int len)
{
    const auto arr_buf = arr.request();
    if (arr_buf.ndim != 1 || arr_buf.shape[0] != len)
        throw_error("All arrays must be one-dimensional with equal size");

    const auto cptr = static_cast<REAL*>(arr_buf.ptr);

    return cptr;
}


// ******************************* Chi parameters ***********************************************

pyArr
chi_photon_wrapper(
    const pyArr& px, const pyArr& py, const pyArr& pz,
    const pyArr& ex, const pyArr& ey, const pyArr& ez,
    const pyArr& bx, const pyArr& by, const pyArr& bz,
    const REAL ref_quantity)
{
    const REAL
        *p_px = nullptr, *p_py = nullptr, *p_pz = nullptr,
        *p_ex = nullptr, *p_ey = nullptr, *p_ez = nullptr,
        *p_bx = nullptr, *p_by = nullptr, *p_bz = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many,
        p_px, p_py, p_pz,
        p_ex, p_ey, p_ez,
        p_bx, p_by, p_bz) =
            check_and_get_pointers(
                px, py, pz,
                ex, ey, ez,
                bx, by, bz);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_phys::chi_photon<REAL, UU>(
                p_px[i], p_py[i], p_pz[i],
                p_ex[i], p_ey[i], p_ez[i],
                p_bx[i], p_by[i], p_bz[i],
                ref_quantity);
    });

    return res;
}

pyArr
chi_ele_pos_wrapper(
    const pyArr& px, const pyArr& py, const pyArr& pz,
    const pyArr& ex, const pyArr& ey, const pyArr& ez,
    const pyArr& bx, const pyArr& by, const pyArr& bz,
    const REAL ref_quantity)
{
    const REAL
        *p_px = nullptr, *p_py = nullptr, *p_pz = nullptr,
        *p_ex = nullptr, *p_ey = nullptr, *p_ez = nullptr,
        *p_bx = nullptr, *p_by = nullptr, *p_bz = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many,
        p_px, p_py, p_pz,
        p_ex, p_ey, p_ez,
        p_bx, p_by, p_bz) =
            check_and_get_pointers(
                px, py, pz,
                ex, ey, ez,
                bx, by, bz);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_phys::chi_ele_pos<REAL, UU>(
                p_px[i], p_py[i], p_pz[i],
                p_ex[i], p_ey[i], p_ez[i],
                p_bx[i], p_by[i], p_bz[i],
                ref_quantity);
    });

    return res;
}

// ______________________________________________________________________________________________



// ******************************* Breit-Wheeler pair production ********************************
using bw_dndt_lookup_table_params =
    pxr_bw::dndt_lookup_table_params<REAL>;

using bw_pair_prod_lookup_table_params =
    pxr_bw::pair_prod_lookup_table_params<REAL>;

using bw_dndt_lookup_table =
    pxr_bw::dndt_lookup_table<REAL, stdVec>;

using bw_pair_prod_lookup_table =
    pxr_bw::pair_prod_lookup_table<REAL, stdVec>;

const auto bw_regular =
    pxr_bw::generation_policy::regular;

const auto bw_force_double =
    pxr_bw::generation_policy::force_internal_double;

pyArr
bw_get_optical_depth_wrapper(
    const pyArr& unf_zero_one_minus_epsi)
{
    const REAL
        *p_unf_zero_one_minus_epsi = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many, p_unf_zero_one_minus_epsi)=
            check_and_get_pointers(unf_zero_one_minus_epsi);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_bw::get_optical_depth<REAL>(
                p_unf_zero_one_minus_epsi[i]);
    });

    return res;
}

pyArr
bw_get_dn_dt_wrapper(
    const pyArr& energy_phot, const pyArr& chi_phot,
    const bw_dndt_lookup_table& ref_table,
    const REAL ref_quantity)
{
    const REAL
        *p_energy_phot = nullptr, *p_chi_phot = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many,
        p_energy_phot, p_chi_phot) =
            check_and_get_pointers(
                energy_phot, chi_phot);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_bw::get_dN_dt<REAL, bw_dndt_lookup_table, UU>(
                p_energy_phot[i], p_chi_phot[i],
                ref_table, ref_quantity);
    });

    return res;
}

void
bw_evolve_optical_depth_wrapper(
    const pyArr& energy_phot, const pyArr& chi_phot,
    const REAL dt, pyArr& optical_depth,
    const bw_dndt_lookup_table& ref_table,
    const REAL ref_quantity)
{
    const REAL
        *p_energy_phot = nullptr, *p_chi_phot = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many,
        p_energy_phot, p_chi_phot) =
            check_and_get_pointers(
                energy_phot, chi_phot);

    auto p_optical_depth =
        check_and_get_pointer_nonconst(optical_depth, how_many);

    PXRQEDPY_FOR(how_many, [&](int i)->void{
        pxr_bw::evolve_optical_depth<REAL, bw_dndt_lookup_table, UU>(
            p_energy_phot[i], p_chi_phot[i],
            dt, p_optical_depth[i],
            ref_table, ref_quantity);
    });
}
// ______________________________________________________________________________________________


// ******************************* Quantum syncrhotron emission *********************************
using qs_dndt_lookup_table_params =
    pxr_qs::dndt_lookup_table_params<REAL>;

using qs_photon_emission_lookup_table_params =
    pxr_qs::photon_emission_lookup_table_params<REAL>;

using qs_dndt_lookup_table =
    pxr_qs::dndt_lookup_table<REAL, stdVec>;

using qs_photon_emission_lookup_table =
    pxr_qs::photon_emission_lookup_table<REAL, stdVec>;

const auto qs_regular =
    pxr_qs::generation_policy::regular;

const auto qs_force_double =
    pxr_qs::generation_policy::force_internal_double;

pyArr
qs_get_optical_depth_wrapper(
    const pyArr& unf_zero_one_minus_epsi)
{
    const REAL
        *p_unf_zero_one_minus_epsi = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many, p_unf_zero_one_minus_epsi)=
            check_and_get_pointers(unf_zero_one_minus_epsi);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_qs::get_optical_depth<REAL>(
                p_unf_zero_one_minus_epsi[i]);
    });

    return res;
}

// ______________________________________________________________________________________________


// ******************************* Schwinger pair production *************************************

pyArr
sc_pair_production_rate_wrapper(
    const pyArr& ex, const pyArr& ey, const pyArr& ez,
    const pyArr& bx, const pyArr& by, const pyArr& bz,
    const REAL ref_quantity)
{
    const REAL
        *p_ex = nullptr, *p_ey = nullptr, *p_ez = nullptr,
        *p_bx = nullptr, *p_by = nullptr, *p_bz = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many,
        p_ex, p_ey, p_ez,
        p_bx, p_by, p_bz) =
            check_and_get_pointers(
                ex, ey, ez,
                bx, by, bz);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_sc::pair_production_rate<REAL, UU>(
                p_ex[i], p_ey[i], p_ez[i],
                p_bx[i], p_by[i], p_bz[i],
                ref_quantity);
    });

    return res;
}

pyArr
sc_expected_pair_number_wrapper(
    const pyArr& ex, const pyArr& ey, const pyArr& ez,
    const pyArr& bx, const pyArr& by, const pyArr& bz,
    const REAL volume, const REAL dt,
    const REAL ref_quantity)
{
    const REAL
        *p_ex = nullptr, *p_ey = nullptr, *p_ez = nullptr,
        *p_bx = nullptr, *p_by = nullptr, *p_bz = nullptr;
    
    size_t how_many = 0;

    std::tie(
        how_many,
        p_ex, p_ey, p_ez,
        p_bx, p_by, p_bz) =
            check_and_get_pointers(
                ex, ey, ez,
                bx, by, bz);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_sc::expected_pair_number<REAL, UU>(
                p_ex[i], p_ey[i], p_ez[i],
                p_bx[i], p_by[i], p_bz[i],
                volume, dt,
                ref_quantity);
    });

    return res;
}

// ______________________________________________________________________________________________


// ******************************* Python module *************************************************


PYBIND11_MODULE(pxr_qed, m) {
    m.doc() = "pybind11 pxr_qed plugin";

    m.attr("PRECISION") = py::str(PXRQEDPY_PRECISION_STRING);
    m.attr("HAS_OPENMP") = py::bool_(PXRQEDPY_OPENMP_FLAG);
    m.attr("UNITS") = py::str(PXRQEDPY_USTRING);

    m.def(
        "chi_photon",
        &chi_photon_wrapper,
        py::arg("px").noconvert(true), py::arg("py").noconvert(true), py::arg("pz").noconvert(true),
        py::arg("ex").noconvert(true), py::arg("ey").noconvert(true), py::arg("ez").noconvert(true),
        py::arg("bx").noconvert(true), py::arg("by").noconvert(true), py::arg("bz").noconvert(true),
        py::arg("ref_quantity") = py::float_(1.0)
        );

    m.def(
        "chi_ele_pos",
        &chi_ele_pos_wrapper,
        py::arg("px").noconvert(true), py::arg("py").noconvert(true), py::arg("pz").noconvert(true),
        py::arg("ex").noconvert(true), py::arg("ey").noconvert(true), py::arg("ez").noconvert(true),
        py::arg("bx").noconvert(true), py::arg("by").noconvert(true), py::arg("bz").noconvert(true),
        py::arg("ref_quantity") = py::float_(1.0)
        );


    auto bw = m.def_submodule( "bw" );

    bw.def(
        "get_optical_depth",
        &bw_get_optical_depth_wrapper,
        "Computes the optical depth of a new photon",
        py::arg("unf_zero_one_minus_epsi").noconvert(true));

    py::class_<bw_dndt_lookup_table_params>(bw,
        "dndt_lookup_table_params",
        "Parameters to generate a dN/dt lookup table")
        .def(py::init<>())
        .def(py::init<REAL,REAL,int>())
        .def("__eq__", &bw_dndt_lookup_table_params::operator==)
        .def_readwrite("chi_phot_min", &bw_dndt_lookup_table_params::chi_phot_min)
        .def_readwrite("chi_phot_max", &bw_dndt_lookup_table_params::chi_phot_max)
        .def_readwrite("chi_phot_how_many", &bw_dndt_lookup_table_params::chi_phot_how_many)
        .def("__repr__",
            [](const bw_dndt_lookup_table_params &a) {
                return 
                    std::string("bw.dndt_lookup_table_params:\n")+
                    std::string("\tchi_phot_min     : ") + float_to_string(a.chi_phot_min)+"\n"+
                    std::string("\tchi_phot_max     : ") + float_to_string(a.chi_phot_max)+"\n"+
                    std::string("\tchi_phot_how_many: ") + std::to_string(a.chi_phot_how_many);
            });

    py::class_<bw_pair_prod_lookup_table_params>(bw,
        "pair_prod_lookup_table_params")
        .def(py::init<>())
        .def(py::init<REAL,REAL,int,int>())
        .def("__eq__", &bw_pair_prod_lookup_table_params::operator==)
        .def_readwrite("chi_phot_min", &bw_pair_prod_lookup_table_params::chi_phot_min)
        .def_readwrite("chi_phot_max", &bw_pair_prod_lookup_table_params::chi_phot_max)
        .def_readwrite("chi_phot_how_many", &bw_pair_prod_lookup_table_params::chi_phot_how_many)
        .def_readwrite("frac_how_many", &bw_pair_prod_lookup_table_params::frac_how_many)
        .def("__repr__",
            [](const bw_pair_prod_lookup_table_params &a) {
                return 
                    std::string("bw.pair_prod_lookup_table_params:\n")+
                    std::string("\tchi_phot_min     : ") + float_to_string(a.chi_phot_min)+"\n"+
                    std::string("\tchi_phot_max     : ") + float_to_string(a.chi_phot_max)+"\n"+
                    std::string("\tchi_phot_how_many: ") + std::to_string(a.chi_phot_how_many)+"\n"+
                    std::string("\tfrac_how_many    : ") + std::to_string(a.frac_how_many);

            });

    py::class_<bw_dndt_lookup_table>(bw,
        "dndt_lookup_table",
        "dN/dt lookup table")
        .def(py::init<>())
        .def(py::init<bw_dndt_lookup_table_params>())
        .def("__eq__", &bw_dndt_lookup_table::operator==)
        .def("generate",
            [&](bw_dndt_lookup_table &self,
                bool do_regular, bool verbose){
                    if(do_regular)
                        self.generate<bw_regular>(verbose);
                    else
                        self.generate<bw_force_double>(verbose);
            },
            py::arg("do_regular") = py::bool_(true),
            py::arg("verbose") = py::bool_(true))
        .def("save_as",
            [&](const bw_dndt_lookup_table &self, const std::string file_name){
                if(!self.is_init())
                    throw_error("Table must be initialized!");
                const auto raw = self.serialize();
                auto of = std::fstream(file_name,
                    std::ios::out | std::ios::binary);
                if( !of )
                    throw_error("Opening file failed!");
                of.write(raw.data(), raw.size());
                of.close();
            },
            py::arg("file_name"))
        .def("load_from",
            [&](bw_dndt_lookup_table &self, const std::string file_name){
                auto input = std::ifstream(file_name,
                    std::ios::ate | std::ios::binary);
                if( !input )
                    throw_error("Opening file failed!");
                const auto pos = input.tellg();
                auto raw = rawVec(pos);
                
                input.seekg(0, std::ios::beg);
                input.read(raw.data(), pos);

                self = bw_dndt_lookup_table{raw};
                input.close();
            },
            py::arg("file_name"))
        .def("interp",
            [&](bw_dndt_lookup_table &self, const pyArr& chi_phot){            
                const REAL* p_chi_phot = nullptr;    
                size_t how_many = 0;
                std::tie(how_many, p_chi_phot)=
                    check_and_get_pointers(chi_phot);

                auto res = pyArr(how_many);
                auto p_res = static_cast<REAL*>(res.request().ptr);

                PXRQEDPY_FOR(how_many, [&](int i){
                    p_res[i] = self.interp(p_chi_phot[i]);
                });
                return res;
            })
        .def("__repr__",
            [](const bw_dndt_lookup_table &a) {
                return 
                    std::string("bw.dndt_lookup_table:\n")+
                    std::string("\tis initialized? : ") + bool_to_string(a.is_init())+"\n";
            });

    py::class_<bw_pair_prod_lookup_table>(bw,
        "pair_prod_lookup_table",
        "Pair production lookup table")
        .def(py::init<>())
        .def(py::init<bw_pair_prod_lookup_table_params>())
        .def("__eq__", &bw_pair_prod_lookup_table::operator==)
        .def("generate",
            [&](bw_pair_prod_lookup_table &self,
                bool do_regular, bool verbose){
                    if(do_regular)
                        self.generate<bw_regular>(verbose);
                    else
                        self.generate<bw_force_double>(verbose);
            },
            py::arg("do_regular") = py::bool_(true),
            py::arg("verbose") = py::bool_(true))
        .def("save_as",
            [&](const bw_pair_prod_lookup_table &self, const std::string file_name){
                if(!self.is_init())
                    throw_error("Table must be initialized!");
                const auto raw = self.serialize();
                auto of = std::fstream(file_name,
                    std::ios::out | std::ios::binary);
                if( !of )
                    throw_error("Opening file failed!");
                of.write(raw.data(), raw.size());
                of.close();
            },
            py::arg("file_name"))
        .def("load_from",
            [&](bw_pair_prod_lookup_table &self, const std::string file_name){
                auto input = std::ifstream(file_name,
                    std::ios::ate | std::ios::binary);
                if( !input )
                    throw_error("Opening file failed!");
                const auto pos = input.tellg();
                auto raw = rawVec(pos);
                
                input.seekg(0, std::ios::beg);
                input.read(raw.data(), pos);

                self = bw_pair_prod_lookup_table{raw};
                input.close();
            },
            py::arg("file_name"))
        .def("interp",
            [&](bw_pair_prod_lookup_table &self,
                    const pyArr& chi_phot, const pyArr& unf_zero_one_minus_epsi){            
                const REAL
                    *p_chi_phot = nullptr, *p_unf_zero_one_minus_epsi = nullptr;
                size_t how_many = 0;
                std::tie(how_many, p_chi_phot, p_unf_zero_one_minus_epsi)=
                    check_and_get_pointers(chi_phot, unf_zero_one_minus_epsi);

                auto res = pyArr(how_many);
                auto p_res = static_cast<REAL*>(res.request().ptr);

                PXRQEDPY_FOR(how_many, [&](int i){
                    p_res[i] = self.interp(p_chi_phot[i], p_unf_zero_one_minus_epsi[i]);
                });
                return res;
            })
        .def("__repr__",
            [](const bw_pair_prod_lookup_table &a) {
                return 
                    std::string("bw.pair_prod_lookup_table:\n")+
                    std::string("\tis initialized? : ") + bool_to_string(a.is_init())+"\n";
            });

    bw.def(
        "get_dn_dt",
        &bw_get_dn_dt_wrapper,
        py::arg("energy_phot").noconvert(true),
        py::arg("chi_phot").noconvert(true),
        py::arg("ref_table"), py::arg("ref_quantity") = py::float_(1.0)
        );

    bw.def(
        "evolve_optical_depth",
        &bw_evolve_optical_depth_wrapper,
        py::arg("energy_phot").noconvert(true),
        py::arg("chi_phot").noconvert(true),
        py::arg("dt"), py::arg("optical_depth").noconvert(true),
        py::arg("ref_table"), py::arg("ref_quantity") = py::float_(1.0)
        );

    auto qs = m.def_submodule( "qs" );

    qs.def(
        "get_optical_depth", 
         &bw_get_optical_depth_wrapper,
        "Computes the optical depth of a new electron or positron",
        py::arg("unf_zero_one_minus_epsi").noconvert(true));

    py::class_<qs_dndt_lookup_table_params>(qs,
        "dndt_lookup_table_params",
        "Parameters to generate a dN/dt lookup table")
        .def(py::init<>())
        .def(py::init<REAL,REAL,int>())
        .def("__eq__", &qs_dndt_lookup_table_params::operator==)
        .def_readwrite("chi_part_min", &qs_dndt_lookup_table_params::chi_part_min)
        .def_readwrite("chi_part_max", &qs_dndt_lookup_table_params::chi_part_max)
        .def_readwrite("chi_part_how_many", &qs_dndt_lookup_table_params::chi_part_how_many)
        .def("__repr__",
            [](const qs_dndt_lookup_table_params &a) {
                return 
                    std::string("qs.dndt_lookup_table_params:\n")+
                    std::string("\tchi_part_min     : ") + float_to_string(a.chi_part_min)+"\n"+
                    std::string("\tchi_part_max     : ") + float_to_string(a.chi_part_max)+"\n"+
                    std::string("\tchi_part_how_many: ") + std::to_string(a.chi_part_how_many);
            });

    py::class_<qs_photon_emission_lookup_table_params>(qs,
        "photon_emission_lookup_table_params")
        .def(py::init<>())
        .def(py::init<REAL,REAL,int,int>())
        .def("__eq__", &qs_photon_emission_lookup_table_params::operator==)
        .def_readwrite("chi_part_min", &qs_photon_emission_lookup_table_params::chi_part_min)
        .def_readwrite("chi_part_max", &qs_photon_emission_lookup_table_params::chi_part_max)
        .def_readwrite("frac_min", &qs_photon_emission_lookup_table_params::frac_min)
        .def_readwrite("chi_part_how_many", &qs_photon_emission_lookup_table_params::chi_part_how_many)
        .def_readwrite("frac_how_many", &qs_photon_emission_lookup_table_params::frac_how_many)
        .def("__repr__",
            [](const qs_photon_emission_lookup_table_params &a) {
                return 
                    std::string("qs.photon_emission_lookup_table_params:\n")+
                    std::string("\tchi_part_min     : ") + float_to_string(a.chi_part_min)+"\n"+
                    std::string("\tchi_part_max     : ") + float_to_string(a.chi_part_max)+"\n"+
                    std::string("\tfrac_min         : ") + float_to_string(a.frac_min)+"\n"+
                    std::string("\tchi_part_how_many: ") + std::to_string(a.chi_part_how_many)+"\n"+
                    std::string("\tfrac_how_many    : ") + std::to_string(a.frac_how_many);
            });

    py::class_<qs_dndt_lookup_table>(qs,
        "dndt_lookup_table",
        "dN/dt lookup table")
        .def(py::init<>())
        .def(py::init<qs_dndt_lookup_table_params>())
        .def("__eq__", &qs_dndt_lookup_table::operator==)
        .def("generate",
            [&](qs_dndt_lookup_table &self,
                bool do_regular, bool verbose){
                    if(do_regular)
                        self.generate<qs_regular>(verbose);
                    else
                        self.generate<qs_force_double>(verbose);
            },
            py::arg("do_regular") = py::bool_(true),
            py::arg("verbose") = py::bool_(true))
        .def("save_as",
            [&](const qs_dndt_lookup_table &self, const std::string file_name){
                if(!self.is_init())
                    throw_error("Table must be initialized!");
                const auto raw = self.serialize();
                auto of = std::fstream(file_name,
                    std::ios::out | std::ios::binary);
                if( !of )
                    throw_error("Opening file failed!");
                of.write(raw.data(), raw.size());
                of.close();
            },
            py::arg("file_name"))
        .def("load_from",
            [&](qs_dndt_lookup_table &self, const std::string file_name){
                auto input = std::ifstream(file_name,
                    std::ios::ate | std::ios::binary);
                if( !input )
                    throw_error("Opening file failed!");
                const auto pos = input.tellg();
                auto raw = rawVec(pos);
                
                input.seekg(0, std::ios::beg);
                input.read(raw.data(), pos);

                self = qs_dndt_lookup_table{raw};
                input.close();
            },
            py::arg("file_name"))
        .def("interp",
            [&](qs_dndt_lookup_table &self, const pyArr& chi_part){            
                const REAL* p_chi_part = nullptr;    
                size_t how_many = 0;
                std::tie(how_many, p_chi_part)=
                    check_and_get_pointers(chi_part);

                auto res = pyArr(how_many);
                auto p_res = static_cast<REAL*>(res.request().ptr);

                PXRQEDPY_FOR(how_many, [&](int i){
                    p_res[i] = self.interp(p_chi_part[i]);
                });
                return res;
            })
        .def("__repr__",
            [](const qs_dndt_lookup_table &a) {
                return 
                    std::string("qs.dndt_lookup_table:\n")+
                    std::string("\tis initialized? : ") + bool_to_string(a.is_init())+"\n";
            });

    py::class_<qs_photon_emission_lookup_table>(qs,
        "photon_emission_lookup_table",
        "Photon emission lookup table")
        .def(py::init<>())
        .def(py::init<qs_photon_emission_lookup_table_params>())
        .def("__eq__", &qs_photon_emission_lookup_table::operator==)
        .def("generate",
            [&](qs_photon_emission_lookup_table &self,
                bool do_regular, bool verbose){
                    if(do_regular)
                        self.generate<qs_regular>(verbose);
                    else
                        self.generate<qs_force_double>(verbose);
            },
            py::arg("do_regular") = py::bool_(true),
            py::arg("verbose") = py::bool_(true))
        .def("save_as",
            [&](const qs_photon_emission_lookup_table &self, const std::string file_name){
                if(!self.is_init())
                    throw_error("Table must be initialized!");
                const auto raw = self.serialize();
                auto of = std::fstream(file_name,
                    std::ios::out | std::ios::binary);
                if( !of )
                    throw_error("Opening file failed!");
                of.write(raw.data(), raw.size());
                of.close();
            },
            py::arg("file_name"))
        .def("load_from",
            [&](qs_photon_emission_lookup_table &self, const std::string file_name){
                auto input = std::ifstream(file_name,
                    std::ios::ate | std::ios::binary);
                if( !input )
                    throw_error("Opening file failed!");
                const auto pos = input.tellg();
                auto raw = rawVec(pos);
                
                input.seekg(0, std::ios::beg);
                input.read(raw.data(), pos);

                self = qs_photon_emission_lookup_table{raw};
                input.close();
            },
            py::arg("file_name"))
        .def("interp",
            [&](qs_photon_emission_lookup_table &self,
                    const pyArr& chi_part, const pyArr& unf_zero_one_minus_epsi){            
                const REAL
                    *p_chi_part = nullptr, *p_unf_zero_one_minus_epsi = nullptr;
                size_t how_many = 0;
                std::tie(how_many, p_chi_part, p_unf_zero_one_minus_epsi)=
                    check_and_get_pointers(chi_part, unf_zero_one_minus_epsi);

                auto res = pyArr(how_many);
                auto p_res = static_cast<REAL*>(res.request().ptr);

                PXRQEDPY_FOR(how_many, [&](int i){
                    p_res[i] = self.interp(p_chi_part[i], p_unf_zero_one_minus_epsi[i]);
                });
                return res;
            })
        .def("__repr__",
            [](const qs_photon_emission_lookup_table &a) {
                return 
                    std::string("qs.pair_prod_lookup_table:\n")+
                    std::string("\tis initialized? : ") + bool_to_string(a.is_init())+"\n";
        });

    auto sc = m.def_submodule( "sc" );

    sc.def("pair_production_rate",
        &sc_pair_production_rate_wrapper,
        "Computes the Schwinger pair production rate using the Nikishov formula",
        py::arg("ex").noconvert(true), py::arg("ey").noconvert(true), py::arg("ez").noconvert(true),
        py::arg("bx").noconvert(true), py::arg("by").noconvert(true), py::arg("bz").noconvert(true),
        py::arg("ref_quantity") = py::float_(1.0)
        );

    sc.def("expected_pair_number",
        &sc_expected_pair_number_wrapper,
        "Computes the expected number of Schwinger pairs using the Nikishov formula",
        py::arg("ex").noconvert(true), py::arg("ey").noconvert(true), py::arg("ez").noconvert(true),
        py::arg("bx").noconvert(true), py::arg("by").noconvert(true), py::arg("bz").noconvert(true),
        py::arg("volume"), py::arg("dt"),
        py::arg("ref_quantity") = py::float_(1.0)
        );


}

// ______________________________________________________________________________________________
