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
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-5;

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double_1 )
{
    const double v = 1.0/3.0;
    const double x = 0.5;
    const double exp = 0.989031074246724;

    const double k_1f3_05 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_1f3_05/exp, double_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_double_2 )
{
    const double v = 1.0/3.0;
    const double x = 1.0;
    const double exp = 0.438430633441534;

    const double k_1f3_1 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_1f3_1/exp , double_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_double_3 )
{
    const double v = 1.0/3.0;
    const double x = 2.0;
    const double exp = 0.116544961296165;

    const double k_1f3_2 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_1f3_2/exp, double_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_double_4 )
{
    const double v = 2.0/3.0;
    const double x = 0.5;
    const double exp = 1.205930464720336;

    const double k_2f3_05 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_2f3_05/exp, double_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_double_5 )
{
    const double v = 2.0/3.0;
    const double x = 1.0;
    const double exp = 0.494475062104208;

    const double k_2f3_1 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_2f3_1/exp, double_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_double_6 )
{
    const double v = 2.0/3.0;
    const double x = 2.0;
    const double exp = 0.124838927488128;

    const double k_2f3_2 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_2f3_2/exp, double_tolerance);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_float_1 )
{
    const float fv = static_cast<float>(1.0/3.0);
    const float fx = static_cast<float>(0.5);
    const float fexp = static_cast<float>(0.989031074246724);

    const float fk_1f3_05 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_1f3_05/fexp, float_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_float_2 )
{
    const float fv = static_cast<float>(1.0/3.0);
    const float fx = static_cast<float>(1.0);
    const float fexp = static_cast<float>(0.438430633441534);

    const float fk_1f3_1 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_1f3_1/fexp, float_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_float_3 )
{
    const float fv = static_cast<float>(1.0/3.0);
    const float fx = static_cast<float>(2.0);
    const float fexp = static_cast<float>(0.116544961296165);

    const float fk_1f3_2 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_1f3_2/fexp, float_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_float_4 )
{
    const float fv = static_cast<float>(2.0/3.0);
    const float fx = static_cast<float>(0.5);
    const float fexp = static_cast<float>(1.205930464720336);

    const float fk_2f3_05 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_2f3_05/fexp, float_tolerance);
}

BOOST_AUTO_TEST_CASE( bessel_functions_float_5 )
{
    const float fv = static_cast<float>(2.0/3.0);
    const float fx = static_cast<float>(1.0);
    const float fexp = static_cast<float>(0.494475062104208);

    const float fk_2f3_1 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_2f3_1/fexp, float_tolerance);
}


BOOST_AUTO_TEST_CASE( bessel_functions_float_6 )
{
    const float fv = static_cast<float>(2.0/3.0);
    const float fx = static_cast<float>(2.0);
    const float fexp = static_cast<float>(0.124838927488128);

    const float fk_2f3_2 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_2f3_2/fexp, float_tolerance);
}
