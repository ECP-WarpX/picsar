#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "picsar_qed/physics/chi_functions.hpp"

#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_core.hpp"
#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"

#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_core.hpp"
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"

#include "picsar_qed/physics/schwinger/schwinger_pair_engine_core.hpp"

#include <vector>
#include <algorithm>
#include <string>
#include <array>
#include <tuple>

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

template <typename Int>
inline int as_int(Int num)
{
    return static_cast<int>(num);
}

using pyArr = py::array_t<REAL>;
using pyBufInfo = py::buffer_info;
using stdVec = std::vector<REAL>;

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

// ______________________________________________________________________________________________


// ******************************* Quantum syncrhotron emission *********************************
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

    auto qs = m.def_submodule( "qs" );

    qs.def(
        "get_optical_depth", 
         &bw_get_optical_depth_wrapper,
        "Computes the optical depth of a new electron or positron",
        py::arg("unf_zero_one_minus_epsi").noconvert(true));


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



#ifdef ABABABABABABABABBA



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
            [](const std::vector<PXRQEDPY_REAL>& unf_zero_one_minus_epsi){
                const auto size = unf_zero_one_minus_epsi.size();
                auto res = std::vector<PXRQEDPY_REAL>(size);          

                PXRQEDPY_FOR(as_int(size), [&](int i){
                    res[i] = pxr_bw::get_optical_depth<PXRQEDPY_REAL>(
                        unf_zero_one_minus_epsi[i]);
                });
            return res;},
        "Computes the optical depth of a new photon");

    bw.def("get_dn_dt", 
            [](const std::vector<PXRQEDPY_REAL>& energy_phot, const std::vector<PXRQEDPY_REAL>& chi_phot,
                const bw_dndt_lookup_table& ref_table, PXRQEDPY_REAL ref_quantity){
                    const auto size = static_cast<int>(chi_phot.size());
                    assert(energy_phot.size() == size);
                    assert(ref_table.is_init());

                    auto res = std::vector<PXRQEDPY_REAL>(size);          
                    
                    PXRQEDPY_FOR(as_int(size), [&](int i){
                         res[i] = pxr_bw::get_dN_dt<
                                PXRQEDPY_REAL, 
                                bw_dndt_lookup_table, 
                                PXRQEDPY_UU>(
                                energy_phot[i], chi_phot[i],
                                ref_table, ref_quantity,nullptr);
                    });
                    return res;                   
                ;}, 
        "Computes dN/dt for Breit-Wheeler pair production");
/*
        bw.def("evolve_optical_depth", 
            [](const std::vector<PXRQEDPY_REAL>& energy_phot, const std::vector<PXRQEDPY_REAL>& chi_phot,
                PXRQEDPY_REAL dt, std::vector<PXRQEDPY_REAL>& optical_depth,
                const bw_dndt_lookup_table& ref_table, PXRQEDPY_REAL ref_quantity){
                    assert(energy_phot.size() == chi_phot.size());
                    assert(optical_depth.size() == chi_phot.size());
                    assert(ref_table.is_init());

                const auto size = static_cast<int>(chi_phot.size());
                PXRQEDPY_FOR(size,[&](int i){

                };),
        "Computes dN/dt for Breit-Wheeler pair production");
        */
// ______________________________________________________________________________________________

}

#endif