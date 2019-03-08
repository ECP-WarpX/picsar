//####### Test module for special functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "special functions"

//Will automatically define a main for this test 
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "special_functions.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Test Bessel functions with double precision
BOOST_AUTO_TEST_CASE( bessel_functions_double )
{
    const double tol = 1.0e-12;
    const double v = 1.0/3.0;
    const double x = 1.0;
    const double exp = 0.438430633441534;

    const double k_1f3_1 = k_v(v,x);

    BOOST_CHECK_SMALL(1.0 - k_1f3_1/exp ,tol );
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( bessel_functions_float )
{
    const float ftol = 1.0e-5;
    const float fv = static_cast<float>(1.0/3.0);
    const float fx = static_cast<float>(1.0);
    const float fexp = static_cast<float>(0.438430633441534);

    const float fk_1f3_1 = k_v(fv,fx);

    BOOST_CHECK_SMALL(1.0f-fk_1f3_1/fexp ,ftol);
}
