//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/core"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include<array>
 #include<utility>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "quantum_sync_engine_core.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-8;
const double double_small = 1e-30;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-3;
const float float_small = 1e-20;

using namespace picsar::multi_physics::phys::quantum_sync;

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

template <typename T>
T constexpr small()
{
    if(std::is_same<T,float>::value)
        return float_small;
    else
        return double_small;
}

template<typename RealType>
struct fake_G_table
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

// ***Test Quantum Synchrotron optical depth assignment
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_core_opt_depth)
{
    check_opt_depth<double>();
    check_opt_depth<float>();
}
