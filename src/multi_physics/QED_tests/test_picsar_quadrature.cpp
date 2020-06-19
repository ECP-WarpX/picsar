//####### Test module for quadrature functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "math/quadrature"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <array>

#include "math_constants.h"
#include "quadrature.hpp"

using namespace picsar::multi_physics::math;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-6;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-3;

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

// ***Test quadrature algorithms in a finite interval

template<typename RealType>
constexpr void test_quadrature_finite_interval()
{
    const auto sin2 = [](RealType x)
        {return sin(x)*sin(x);};

    const auto a = static_cast<RealType>(0.0);
    const auto b = static_cast<RealType>(2.0)*pi<RealType>;

    const auto exp_res = pi<RealType>;

    const auto res_a_b = quad_a_b<RealType>(sin2, a, b);
    const auto res_a_b_s = quad_a_b_s<RealType>(sin2, a, b);
    const auto res_tanh_sinh =
        generic_quad_a_b<RealType, quadrature_algorithm::tanh_sinh>(
            sin2, a, b);
    const auto res_gauss_kronrod15 =
        generic_quad_a_b<RealType, quadrature_algorithm::gauss_kronrod15>(
            sin2, a, b);
    const auto res_gauss_kronrod31 =
        generic_quad_a_b<RealType, quadrature_algorithm::gauss_kronrod31>(
            sin2, a, b);
    const auto res_gauss_kronrod41 =
        generic_quad_a_b<RealType, quadrature_algorithm::gauss_kronrod41>(
            sin2, a, b);
    const auto res_gauss_kronrod51 =
        generic_quad_a_b<RealType, quadrature_algorithm::gauss_kronrod51>(
            sin2, a, b);
    const auto res_gauss_kronrod61 =
        generic_quad_a_b<RealType, quadrature_algorithm::gauss_kronrod61>(
            sin2, a, b);

    BOOST_CHECK_SMALL((res_a_b-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_a_b_s-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_tanh_sinh-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_gauss_kronrod15-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_gauss_kronrod31-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_gauss_kronrod41-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_gauss_kronrod51-exp_res)/exp_res, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_gauss_kronrod61-exp_res)/exp_res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_quadrature_finite_interval )
{
    test_quadrature_finite_interval <double>();
    test_quadrature_finite_interval <float>();
}

// *******************************

// ***Test quadrature algorithms in a semi-infinite interval

template<typename RealType>
constexpr void test_quadrature_infinite_interval_1()
{
    const auto expx2 = [](RealType x)
        {return exp(-x*x);};

    const auto a = static_cast<RealType>(0.0);

    const auto exp_res = static_cast<RealType>(sqrt(pi<double>)/2.0);

    const auto res_expx2 = quad_a_inf<RealType>(expx2, a);

    BOOST_CHECK_SMALL((res_expx2-exp_res)/exp_res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_quadrature_infinite_interval_1 )
{
    test_quadrature_infinite_interval_1<double>();
    test_quadrature_infinite_interval_1<float>();
}


template<typename RealType>
constexpr void test_quadrature_infinite_interval_2()
{
    const auto datan = [](RealType x)
        {return 1.0/(1.0+x*x);};

    const auto a = static_cast<RealType>(0.0);

    const auto exp_res = static_cast<RealType>(0.5)*pi<RealType>;

    const auto res_datan = quad_a_inf<RealType>(datan, a);

    BOOST_CHECK_SMALL((res_datan-exp_res)/exp_res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_quadrature_infinite_interval_2 )
{
    test_quadrature_infinite_interval_2<double>();
    test_quadrature_infinite_interval_2<float>();
}

// *******************************
