//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "vector functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 //Include Boost unit tests library & library for floating point comparison
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>

#include "vec_functions.hpp"

using namespace picsar::multi_physics;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-13;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-5;

// ------------- Tests --------------

//Test norm in double precision
BOOST_AUTO_TEST_CASE( vec_functions_norm_double_1 )
{
    vec3<double> vv{1.0,-2.0,3.0};
    double exp = sqrt(1.0 + 4.0 + 9.0);
    BOOST_CHECK_SMALL(1.0 - norm(vv)/exp, double_tolerance);
}

//Test norm in single precision
BOOST_AUTO_TEST_CASE( vec_functions_norm_single_1 )
{
    vec3<float> fvv{1.0f,-2.0f,3.0f};
    float fexp = sqrt(1.0f + 4.0f + 9.0f);
    BOOST_CHECK_SMALL(1.0f - norm(fvv)/fexp, float_tolerance);
}

//Test dot in double precision
BOOST_AUTO_TEST_CASE( vec_functions_dot_double_1 )
{
    vec3<double> vv1{1.0,-2.0,3.0};
    vec3<double> vv2{-1.0, 0.0,5.0};
    double exp = 14.0;
    BOOST_CHECK_SMALL(1.0 - dot(vv1, vv2)/exp, double_tolerance);
}

//Test dot in single precision
BOOST_AUTO_TEST_CASE( vec_functions_dot_single_1 )
{
    vec3<float> vv1{1.0f,-2.0f,3.0f};
    vec3<float> vv2{-1.0f, 0.0f,5.0f};
    float exp = 14.0f;
    BOOST_CHECK_SMALL(1.0f - dot(vv1, vv2)/exp, float_tolerance);
}

//Test cross in double precision
BOOST_AUTO_TEST_CASE( vec_functions_cross_double_1 )
{
    vec3<double> vv1{1./3., -1./4., 1./5.};
    vec3<double> vv2{-2./3., -3./4., 4./5.};
    vec3<double> exp{-1./20., -2./5., -5./12.};

    vec3<double> res = cross(vv1,vv2);
    BOOST_CHECK_SMALL(1.0 - res[0]/exp[0], double_tolerance);
    BOOST_CHECK_SMALL(1.0 - res[1]/exp[1], double_tolerance);
    BOOST_CHECK_SMALL(1.0 - res[2]/exp[2], double_tolerance);
}

//Test cross in single precision
BOOST_AUTO_TEST_CASE( vec_functions_cross_single_1 )
{
    vec3<float> fvv1{1.0f/3.0f, -1.0f/4.0f, 1.0f/5.0f};
    vec3<float> fvv2{-2.0f/3.0f, -3.0f/4.0f, 4.0f/5.0f};
    vec3<float> fexp{-1.0f/20.0f, -2.0f/5.0f, -5.0f/12.};

    vec3<float> fres = cross(fvv1,fvv2);
    BOOST_CHECK_SMALL(1.0f - fres[0]/fexp[0], float_tolerance);
    BOOST_CHECK_SMALL(1.0f - fres[1]/fexp[1], float_tolerance);
    BOOST_CHECK_SMALL(1.0f - fres[2]/fexp[2], float_tolerance);
}
