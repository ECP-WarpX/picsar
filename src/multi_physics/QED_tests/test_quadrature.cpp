//####### Test module for quadrature ####################################

//Define Module name
 #define BOOST_TEST_MODULE "quadratures"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 //Include Boost unit tests library & library for floating point comparison
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>

#include "quadrature.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-6;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-3;

//Test quadrature with double precision
BOOST_AUTO_TEST_CASE( quadrature_finite_interval_double_1 )
{
    auto sin2 = [](double x)
        {return sin(x)*sin(x);};

    const double a = 0.0;
    const double b = 2.0*pi;

    const double exp_res_sin2 = pi;

    double res_sin2 = quad_a_b<double>(sin2, a, b);

    BOOST_CHECK_SMALL(1.0 - res_sin2/exp_res_sin2, double_tolerance);
}

BOOST_AUTO_TEST_CASE( quadrature_infinite_interval_double_1 )
{
    auto expx2 = [](double x)
        {return exp(-x*x);};

    const double a = 0.0;

    const double exp_expx2 = sqrt(pi)/2.0;

    double res_expx2 = quad_a_inf<double>(expx2, a);

    BOOST_CHECK_SMALL(1.0 - res_expx2/exp_expx2, double_tolerance);
}

BOOST_AUTO_TEST_CASE( quadrature_infinite_interval_double_2 )
{
    auto datan = [](double x)
        {return 1.0/(1.0+x*x);};

    const double a = 0.0;

    const double exp_datan = pi/2.0;

    double res_datan = quad_a_inf<double>(datan, a);

    BOOST_CHECK_SMALL(1.0 - res_datan/exp_datan, double_tolerance);
}

//Test Bessel functions with single precision
BOOST_AUTO_TEST_CASE( quadrature_finite_interval_float_1 )
{
    auto fsin2 = [](float x)
        {return sin(x)*sin(x);};

        const float fa = static_cast<float>(0.0);
        const float fb = static_cast<float>(2.0*pi);

        const float fexp_res_sin2 = static_cast<float>(pi);

        float fres_sin2 = quad_a_b<float>(fsin2, fa, fb);

        BOOST_CHECK_SMALL(1.0f - fres_sin2/fexp_res_sin2, float_tolerance);
}

BOOST_AUTO_TEST_CASE( quadrature_infinite_interval_float_1 )
{
    auto fexpx2 = [](float x)
        {return exp(-x*x);};

    const float fa = static_cast<float>(0.0);

    const float fexp_expx2 = static_cast<float>(sqrt(pi)/2.0);

    float fres_expx2 = quad_a_inf<float>(fexpx2, fa);

    BOOST_CHECK_SMALL(1.0f - fres_expx2/fexp_expx2, float_tolerance);
}

BOOST_AUTO_TEST_CASE( quadrature_infinite_interval_float_2 )
{
    auto fdatan = [](float x)
        {return 1.0f/(1.0f+x*x);};

    const float fa = static_cast<float>(0.0);

    const float fexp_datan = static_cast<float>(pi/2.0);

    float fres_datan = quad_a_inf<float>(fdatan, fa);

    BOOST_CHECK_SMALL(1.0f - fres_datan/fexp_datan, float_tolerance);
}
