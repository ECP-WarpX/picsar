//####### Test module for special functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "special functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS

#include "special_functions.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-10;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//Templated tolerance
template <typename T>
T tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//Test Bessel functions generic
template<typename T, typename WHATEVER>
void bessel_functions_test(WHATEVER _v, WHATEVER _x, WHATEVER _exp)
{
    const T v = static_cast<T>(_v);
    const T x = static_cast<T>(_x);
    const T exp = static_cast<T>(_exp);

    T res = k_v(v,x);
    BOOST_CHECK_SMALL((res-exp)/exp, tolerance<T>());
}


//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_1 )
{
    bessel_functions_test<double>(1.0/3.0, 0.5, 0.989031074246724);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_single_1 )
{
    bessel_functions_test<float>(1.0/3.0, 0.5, 0.989031074246724);
}

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_2 )
{
    bessel_functions_test<double>(1.0/3.0, 1.0, 0.438430633441534);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_single_2 )
{
    bessel_functions_test<float>(1.0/3.0, 1.0, 0.438430633441534);
}

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_3 )
{
    bessel_functions_test<double>(1.0/3.0, 2.0, 0.116544961296165);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_single_3 )
{
    bessel_functions_test<float>(1.0/3.0, 2.0, 0.116544961296165);
}

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_4 )
{
    bessel_functions_test<double>(2.0/3.0, 0.5, 1.205930464720336);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_float_4 )
{
    bessel_functions_test<float>(2.0/3.0, 0.5, 1.205930464720336);
}

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_5 )
{
    bessel_functions_test<double>(2.0/3.0, 1.0, 0.494475062104208);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_float_5 )
{
    bessel_functions_test<float>(2.0/3.0, 1.0, 0.494475062104208);
}

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_6 )
{
    bessel_functions_test<double>(2.0/3.0, 2.0, 0.124838927488128);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_float_6 )
{
    bessel_functions_test<float>(2.0/3.0, 2.0, 0.124838927488128);
}
