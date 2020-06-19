//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/breit_wheeler/core"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include<array>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "breit_wheeler_engine_core.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-8;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-3;

using namespace picsar::multi_physics::phys::breit_wheeler;

using namespace picsar::multi_physics::phys;

using namespace picsar::multi_physics::math;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

template<typename RealType>
struct fake_T_table
{
    RealType interp(RealType chi) const {
        m_chi = chi;
        return static_cast<RealType>(m_res);
    }

    RealType m_res;
    mutable RealType m_chi;
};

// ------------- Tests --------------
template <typename RealType>
void check_opt_depth()
{
    const auto fake_rand_zero_one_minus_epsi =
        std::array<double, 11>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999};

    for (auto ff : fake_rand_zero_one_minus_epsi){
        RealType res = get_optical_depth(ff);
        RealType expected = -log(static_cast<RealType>(1.0) - ff);
        BOOST_CHECK_EQUAL(res, expected);
    }

}

// ***Test Breit Wheeler dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_core_opt_depth)
{
    check_opt_depth<double>();
    check_opt_depth<float>();
}

template <typename RealType, unit_system UnitSystem>
void check_dndt(RealType ref_q = one<RealType>)
{
    fake_T_table<RealType> fake_table;

    const auto Tvals = std::array<RealType,4> {1e-6,1e-4,2e-1};

    BOOST_CHECK_EQUAL(0,0);
}

// ***Test Breit Wheeler dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_core_get_dndt)
{
    const double reference_length = 800.0e-9;
    const double reference_omega = 2.0*pi<double>*light_speed<double>/
        reference_length;

    check_dndt<double, unit_system::SI>();
    check_dndt<double, unit_system::norm_omega>(reference_omega);
    check_dndt<double, unit_system::norm_lambda>(reference_length);
    check_dndt<double, unit_system::heaviside_lorentz>();
    check_dndt<float, unit_system::SI>();
    check_dndt<float, unit_system::norm_omega>(reference_omega);
    check_dndt<float, unit_system::norm_lambda>(reference_length);
    check_dndt<float, unit_system::heaviside_lorentz>();
}