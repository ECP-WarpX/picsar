//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "math/special_functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <array>

#include "special_functions.hpp"

using namespace picsar::multi_physics::math;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

// ------------- Tests --------------

//***Test Bessel functions

template<typename T, typename WHATEVER>
constexpr void bessel_functions_test(WHATEVER t_v, WHATEVER t_x, WHATEVER t_exp)
{
    const T v = static_cast<T>(t_v);
    const T x = static_cast<T>(t_x);
    const T exp = static_cast<T>(t_exp);

    T res = k_v(v,x);
    BOOST_CHECK_SMALL((res-exp)/exp, tolerance<T>());
}

template<typename T>
constexpr void test_case()
{
    bessel_functions_test<T>(1.0/3.0, 0.5, 0.989031074246724);
    bessel_functions_test<T>(1.0/3.0, 1.0, 0.438430633441534);
    bessel_functions_test<T>(1.0/3.0, 2.0, 0.116544961296165);
    bessel_functions_test<T>(2.0/3.0, 0.5, 1.205930464720336);
    bessel_functions_test<T>(2.0/3.0, 1.0, 0.494475062104208);
    bessel_functions_test<T>(2.0/3.0, 2.0, 0.124838927488128);
}

BOOST_AUTO_TEST_CASE( picsar_bessel_functions )
{
    test_case<double>();
    test_case<float>();
}

//*******************************
