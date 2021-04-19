/**
* This file contains the python bindings for the PICSAR QED library
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#ifdef PXQEDPY_HAS_OPENMP
    #include <omp.h>
#endif

#include "picsar_qed/math/vec_functions.hpp"

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
namespace pxr_math = picsar::multi_physics::math;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;
namespace pxr_sc = picsar::multi_physics::phys::schwinger;

#ifdef PXRQEDPY_DOUBLE_PRECISION
    using REAL = double;
    const auto PXRQEDPY_PRECISION_STRING = std::string{"double"};
#else
    using REAL = float;
    const auto PXRQEDPY_PRECISION_STRING = std::string{"single"};
#endif

#ifdef PXQEDPY_HAS_OPENMP
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

#if defined(PXRQEDPY_SI)
    const auto UU = pxr_phys::unit_system::SI;
    const auto PXRQEDPY_USTRING = std::string{"SI"};
#elif defined(PXRQEDPY_NORM_OMEGA)
    const auto UU = pxr_phys::unit_system::norm_omega;
    const auto PXRQEDPY_USTRING = std::string{"NORM_OMEGA"};
#elif defined(PXRQEDPY_NORM_LAMBDA)
    const auto UU = pxr_phys::unit_system::norm_lambda;
    const auto PXRQEDPY_USTRING = std::string{"NORM_LAMBDA"};
#elif defined(PXRQEDPY_UNITS)
    const auto UU = pxr_phys::unit_system::heaviside_lorentz;
    const auto PXRQEDPY_USTRING = std::string{"HEAVISIDE_LORENTZ"};
#else
    #error Incorrect units choice!
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

    PXRQEDPY_FOR(how_many, [&](int i){
        pxr_bw::evolve_optical_depth<REAL, bw_dndt_lookup_table, UU>(
            p_energy_phot[i], p_chi_phot[i],
            dt, p_optical_depth[i],
            ref_table, ref_quantity);
    });
}

auto
bw_generate_breit_wheeler_pairs_wrapper(
    const pyArr& chi_phot,
    const pyArr& phot_px, const pyArr& phot_py, const pyArr& phot_pz,
    const pyArr& unf_zero_one_minus_epsi,
    const bw_pair_prod_lookup_table& ref_table,
    const REAL ref_quantity)
{
    const REAL
        *p_chi_phot = nullptr,
        *p_phot_px = nullptr, *p_phot_py = nullptr, *p_phot_pz = nullptr,
        *p_unf_zero_one_minus_epsi;

    size_t how_many = 0;

    std::tie(
        how_many,
        p_chi_phot, p_phot_px, p_phot_py, p_phot_pz,
        p_unf_zero_one_minus_epsi) =
            check_and_get_pointers(
                chi_phot, phot_px, phot_py, phot_pz, unf_zero_one_minus_epsi);

    auto ele_px = pyArr(how_many);
    auto ele_py = pyArr(how_many);
    auto ele_pz = pyArr(how_many);
    auto pos_px = pyArr(how_many);
    auto pos_py = pyArr(how_many);
    auto pos_pz = pyArr(how_many);
    auto p_ele_px = static_cast<REAL*>(ele_px.request().ptr);
    auto p_ele_py = static_cast<REAL*>(ele_py.request().ptr);
    auto p_ele_pz = static_cast<REAL*>(ele_pz.request().ptr);
    auto p_pos_px = static_cast<REAL*>(pos_px.request().ptr);
    auto p_pos_py = static_cast<REAL*>(pos_py.request().ptr);
    auto p_pos_pz = static_cast<REAL*>(pos_pz.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        auto ele_mom = pxr_math::vec3<REAL>{};
        auto pos_mom = pxr_math::vec3<REAL>{};
        pxr_bw::generate_breit_wheeler_pairs<REAL, bw_pair_prod_lookup_table, UU>(
            p_chi_phot[i],
            pxr_math::vec3<REAL>{p_phot_px[i], p_phot_py[i], p_phot_pz[i]},
            p_unf_zero_one_minus_epsi[i],
            ref_table,
            ele_mom, pos_mom,
            ref_quantity);

        p_ele_px[i] = ele_mom[0];
        p_ele_py[i] = ele_mom[1];
        p_ele_pz[i] = ele_mom[2];

        p_pos_px[i] = pos_mom[0];
        p_pos_py[i] = pos_mom[1];
        p_pos_pz[i] = pos_mom[2];
    });

    return std::make_tuple(
        std::move(ele_px),std::move(ele_py), std::move(ele_pz),
        std::move(pos_px),std::move(pos_py), std::move(pos_pz));
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

pyArr
qs_get_dn_dt_wrapper(
    const pyArr& energy_part, const pyArr& chi_part,
    const qs_dndt_lookup_table& ref_table,
    const REAL ref_quantity)
{
    const REAL
        *p_energy_part = nullptr, *p_chi_part = nullptr;

    size_t how_many = 0;

    std::tie(
        how_many,
        p_energy_part, p_chi_part) =
            check_and_get_pointers(
                energy_part, chi_part);

    auto res = pyArr(how_many);
    auto p_res = static_cast<REAL*>(res.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        p_res[i] =
            pxr_qs::get_dN_dt<REAL, qs_dndt_lookup_table, UU>(
                p_energy_part[i], p_chi_part[i],
                ref_table, ref_quantity);
    });

    return res;
}

void
qs_evolve_optical_depth_wrapper(
    const pyArr& energy_part, const pyArr& chi_part,
    const REAL dt, pyArr& optical_depth,
    const qs_dndt_lookup_table& ref_table,
    const REAL ref_quantity)
{
    const REAL
        *p_energy_part = nullptr, *p_chi_part = nullptr;

    size_t how_many = 0;

    std::tie(
        how_many,
        p_energy_part, p_chi_part) =
            check_and_get_pointers(
                energy_part, chi_part);

    auto p_optical_depth =
        check_and_get_pointer_nonconst(optical_depth, how_many);

    PXRQEDPY_FOR(how_many, [&](int i){
        pxr_qs::evolve_optical_depth<REAL, qs_dndt_lookup_table, UU>(
            p_energy_part[i], p_chi_part[i],
            dt, p_optical_depth[i],
            ref_table, ref_quantity);
    });
}

auto
qs_generate_photon_update_momentum_wrapper(
    const pyArr& chi_part,
    pyArr& part_px, pyArr& part_py, pyArr& part_pz,
    const pyArr& unf_zero_one_minus_epsi,
    const qs_photon_emission_lookup_table& ref_table,
    const REAL ref_quantity)
{
    const REAL
        *p_chi_part = nullptr, *p_unf_zero_one_minus_epsi = nullptr;

    size_t how_many = 0;

    std::tie(
        how_many,
        p_chi_part, p_unf_zero_one_minus_epsi) =
            check_and_get_pointers(chi_part, unf_zero_one_minus_epsi);

    auto p_part_px =
        check_and_get_pointer_nonconst(part_px, how_many);
    auto p_part_py =
        check_and_get_pointer_nonconst(part_py, how_many);
    auto p_part_pz =
        check_and_get_pointer_nonconst(part_pz, how_many);

    auto phot_px = pyArr(how_many);
    auto phot_py = pyArr(how_many);
    auto phot_pz = pyArr(how_many);
    auto p_phot_px = static_cast<REAL*>(phot_px.request().ptr);
    auto p_phot_py = static_cast<REAL*>(phot_py.request().ptr);
    auto p_phot_pz = static_cast<REAL*>(phot_pz.request().ptr);

    PXRQEDPY_FOR(how_many, [&](int i){
        auto part_mom = pxr_math::vec3<REAL>{p_part_px[i], p_part_py[i], p_part_pz[i]};
        auto phot_mom = pxr_math::vec3<REAL>{};
        pxr_qs::generate_photon_update_momentum<REAL, qs_photon_emission_lookup_table, UU>(
            p_chi_part[i],
            part_mom,
            p_unf_zero_one_minus_epsi[i],
            ref_table,
            phot_mom,
            ref_quantity);

        p_part_px[i] = part_mom[0];
        p_part_py[i] = part_mom[1];
        p_part_pz[i] = part_mom[2];
        p_phot_px[i] = phot_mom[0];
        p_phot_py[i] = phot_mom[1];
        p_phot_pz[i] = phot_mom[2];
    });

    return std::make_tuple(
        std::move(phot_px),std::move(phot_py), std::move(phot_pz));
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

    bw.def(
        "generate_breit_wheeler_pairs",
        &bw_generate_breit_wheeler_pairs_wrapper,
        py::arg("chi_phot").noconvert(true),
        py::arg("phot_px").noconvert(true),
        py::arg("phot_py").noconvert(true),
        py::arg("phot_pz").noconvert(true),
        py::arg("unf_zero_one_minus_epsi").noconvert(true),
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

    qs.def(
        "get_dn_dt",
        &qs_get_dn_dt_wrapper,
        py::arg("energy_part").noconvert(true),
        py::arg("chi_part").noconvert(true),
        py::arg("ref_table"), py::arg("ref_quantity") = py::float_(1.0)
        );

    qs.def(
        "evolve_optical_depth",
        &qs_evolve_optical_depth_wrapper,
        py::arg("energy_part").noconvert(true),
        py::arg("chi_part").noconvert(true),
        py::arg("dt"), py::arg("optical_depth").noconvert(true),
        py::arg("ref_table"), py::arg("ref_quantity") = py::float_(1.0)
        );

    qs.def(
        "generate_photon_update_momentum",
        &qs_generate_photon_update_momentum_wrapper,
        py::arg("chi_part").noconvert(true),
        py::arg("part_px").noconvert(true),
        py::arg("part_py").noconvert(true),
        py::arg("part_pz").noconvert(true),
        py::arg("unf_zero_one_minus_epsi").noconvert(true),
        py::arg("ref_table"), py::arg("ref_quantity") = py::float_(1.0)
        );

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
