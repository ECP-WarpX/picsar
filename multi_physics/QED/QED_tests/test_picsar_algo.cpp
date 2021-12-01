//####### Test module for picsar_algo ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/picsar_algo"

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

//Force the use of PICSAR implementation of upper_bound and lower_bound for debug purposes
#define PXRMP_PICSAR_UPPER_BOUND
#define PXRMP_PICSAR_LOWER_BOUND
#include <picsar_qed/utils/picsar_algo.hpp>

#include <array>
#include <algorithm>

using namespace picsar::multi_physics::utils;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-6;

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

// ***Test upper_bound

BOOST_AUTO_TEST_CASE( picsar_upper_bound_1 )
{
    const auto arr = std::array<double,5>{0.0, 1.0, 2.0, 3.0, 4.0};

    for (const auto nn : std::vector<double>{-1.,0.,0.1,1,1.1,3.9,4.0,5.0}){
        BOOST_CHECK_EQUAL(
            picsar_upper_bound(arr.begin(), arr.end(), nn),
            std::upper_bound(arr.begin(), arr.end(), nn));
    }
}

// ***Test lower_bound

BOOST_AUTO_TEST_CASE( picsar_lower_bound_1 )
{
    const auto arr = std::array<double,5>{0.0, 1.0, 2.0, 3.0, 4.0};

    for (const auto nn : std::vector<double>{-1.,0.,0.1,1,1.1,3.9,4.0,5.0}){
        BOOST_CHECK_EQUAL(
            picsar_lower_bound(arr.begin(), arr.end(), nn),
            std::lower_bound(arr.begin(), arr.end(), nn));
    }
}

// ***Test upper_bound_functor

BOOST_AUTO_TEST_CASE( picsar_upper_bound_functor_1 )
{
    const auto arr = std::array<double,5>{0.0, 1.0, 2.0, 3.0, 4.0};

    for (const auto nn : std::vector<double>{-1.,0.,0.1,1,1.1,3.9,4.0,5.0}){
        BOOST_CHECK_EQUAL(
            picsar_upper_bound_functor(0, arr.size(), nn, [&](size_t i){
                return arr[i];}),
            std::distance(
                arr.begin(),
                std::upper_bound(arr.begin(), arr.end(), nn)));

    }
}

// ***Test lower_bound_functor

BOOST_AUTO_TEST_CASE( picsar_lower_bound_functor_1 )
{
    const auto arr = std::array<double,5>{0.0, 1.0, 2.0, 3.0, 4.0};

    for (const auto nn : std::vector<double>{-1.,0.,0.1,1,1.1,3.9,4.0,5.0}){
        BOOST_CHECK_EQUAL(
            picsar_lower_bound_functor(0, arr.size(), nn, [&](size_t i){
                return arr[i];}),
            std::distance(
                arr.begin(),
                std::lower_bound(arr.begin(), arr.end(), nn)));

    }
}


// ***Test linear_interp

template<typename RealType>
void check_linear_interp()
{
    const auto linear_func = [](RealType x){
        return static_cast<RealType>(11.3 - 7.9*x);
    };

    const RealType x0 = 3.1;
    const RealType x1 = 4.9;
    const RealType y0 = linear_func(x0);
    const RealType y1 = linear_func(x1);
    const RealType x = 4.1;
    const RealType yexp = linear_func(x);
    const RealType y = linear_interp(x0, x1, y0, y1, x);

    BOOST_CHECK_SMALL((y-yexp)/y, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_linear_interp )
{
    check_linear_interp<double>();
    check_linear_interp<float>();
}

// *******************************

// ***Test bilinear_interp

template<typename RealType>
void check_bilinear_interp()
{
    const auto bilinear_func = [](RealType x, RealType y){
        return static_cast<RealType>(11.3 - 7.9*x + 5.0*y);
    };

    const RealType x0 = 3.1;
    const RealType x1 = 4.9;
    const RealType y0 = -1.9;
    const RealType y1 = -1.0;
    const RealType z00 = bilinear_func(x0,y0);
    const RealType z01 = bilinear_func(x0,y1);
    const RealType z10 = bilinear_func(x1,y0);
    const RealType z11 = bilinear_func(x1,y1);
    const RealType x = 4.1;
    const RealType y = -1.5;
    const RealType zexp = bilinear_func(x,y);
    const RealType z = bilinear_interp(
        x0, x1, y0, y1, z00, z01, z10, z11, x, y);

    BOOST_CHECK_SMALL((z-zexp)/z, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_bilinear_interp )
{
    check_bilinear_interp<double>();
    check_bilinear_interp<float>();
}

// *******************************
