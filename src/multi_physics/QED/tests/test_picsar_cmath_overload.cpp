//####### Test module for the overload of cmath functions ##################

//Define Module name
 #define BOOST_TEST_MODULE "math/cmath_overload"

#include <cmath>
#include <vector>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#define PXRMP_PREVENT_USE_STD_FOR_MATH
#include "cmath_overloads.hpp"

using namespace picsar::multi_physics::math;

// ------------- Tests --------------

// ***Test cmath overloads

template<typename RealType>
void test_case_cmath_overloads()
{
    const auto vals = std::vector<RealType>{1.0e-30, 1.0e-20, 1.0e-10, 1.0,
     3.141459, 1.0e10, 1.0e20, 1.0e30};

     for (const auto val : vals){
        BOOST_CHECK_EQUAL(m_sqrt(val), std::sqrt(val));
        BOOST_CHECK_EQUAL(m_cbrt(val), std::cbrt(val));
        BOOST_CHECK_EQUAL(m_exp(val), std::exp(val));
        BOOST_CHECK_EQUAL(m_exp(-val), std::exp(-val));
        BOOST_CHECK_EQUAL(m_log(val), std::log(val));
        BOOST_CHECK_EQUAL(m_tanh(val), std::tanh(val));
        BOOST_CHECK_EQUAL(m_tanh(-val), std::tanh(-val));
     }

    const auto vals_floor_fabs = std::vector<RealType>{0.0, 0.001, 0.3, 0.7,
         1.11, 1.5, 1.55, 20.9, 100.56, 1000.24};

    for (const auto val : vals_floor_fabs){
        BOOST_CHECK_EQUAL(m_floor(val), std::floor(val));
        BOOST_CHECK_EQUAL(m_floor(-val), std::floor(-val));

        BOOST_CHECK_EQUAL(m_fabs(val), std::fabs(val));
        BOOST_CHECK_EQUAL(m_fabs(-val), std::fabs(-val));
    }
}

BOOST_AUTO_TEST_CASE( picsar_cmath_overloads )
{
    test_case_cmath_overloads<double>();
    test_case_cmath_overloads<float>();
}

// *******************************
