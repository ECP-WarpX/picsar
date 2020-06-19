//####### Test module for schwinger engine (core functions) ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/schwinger"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <vector>
#include <array>

#include "schwinger_pair_engine_core.hpp"

using namespace picsar::multi_physics::utils;

using namespace picsar::multi_physics::phys;

using namespace picsar::multi_physics::math;

using namespace picsar::multi_physics::phys::schwinger;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-8;

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

std::vector<std::array<double,3>> E =
    {std::array<double,3>{-3.96271330e+19,  3.23971426e+19, -3.64982966e+19},
    std::array<double,3>{ 1.56378273e+18 ,-2.75598844e+18 ,-6.36294219e+18},
    std::array<double,3>{-4.39048913e+17 ,-2.51103888e+17 ,-2.76438872e+17},
    std::array<double,3>{2.29889814e+16 , 3.61624665e+15 , 4.44570016e+16},
    std::array<double,3>{ 4.75761748e+15 , 2.94590967e+15, -4.36371026e+15},
    std::array<double,3>{-2.55980495e+17 ,-5.14196104e+16 , 6.25509700e+17},
    std::array<double,3>{0. , 0. ,  0.}
};

std::vector<std::array<double,3>> B =
    {std::array<double,3>{1.83642759e+10 , 1.34563931e+11, -3.16337166e+10},
    std::array<double,3>{ 5.94495226e+09 ,-2.37569313e+09 ,-6.25682892e+09},
    std::array<double,3>{-1.99973282e+09 ,-1.45199067e+09 ,-1.75951091e+09},
    std::array<double,3>{ 8.50494135e+07 , 1.29804272e+08 ,-1.20161632e+08},
    std::array<double,3>{ 4722566.13842951,  7784548.86685482, -7840892.65626958},
    std::array<double,3>{ 0.0 , 0.0, 0.0},
    std::array<double,3>{8.09766737e+08 , 8.47995019e+08 , 2.04094773e+09}
};

std::vector<double> res_exp =
    {2.8362320807767248e+17,
    2276477854308499.5,
    70319852141.39209,
    1.568366399483277e-68,
    6.270694104358953e-249,
    61874419984.33997,
    0.0};

double volume = (1.0e-27);
double dt = 1.0e-15;

template<unit_system UnitSystem, typename RealType>
void test_expected_pair_number(RealType ref = zero<RealType>)
{
    for(int i = 0 ; i < E.size(); ++i){
        const auto fe = conv<quantity::E, unit_system::SI, UnitSystem, RealType>::
            fact(1.0, ref);
        const auto fb = conv<quantity::B, unit_system::SI, UnitSystem, RealType>::
            fact(1.0, ref);
        const auto fv = conv<quantity::volume, unit_system::SI, UnitSystem, RealType>::
            fact(1.0, ref);
        const auto ft = conv<quantity::time, unit_system::SI, UnitSystem, RealType>::
            fact(1.0, ref);
        const auto res = expected_pair_number<RealType, UnitSystem>(
            E[i][0]*fe, E[i][1]*fe, E[i][2]*fe, B[i][0]*fb, B[i][1]*fb, B[i][2]*fb,
            volume*fv, dt*ft, ref);

        if(res_exp[i] <= tolerance<RealType>())
            BOOST_CHECK_SMALL(res, tolerance<RealType>());
        else
            BOOST_CHECK_SMALL((res-static_cast<RealType>(res_exp[i]))/static_cast<RealType>(res_exp[i]), tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_schwinger_core_expected_pair_number )
{
    const double reference_length = 800.0e-9;
    const double reference_omega = 2.0*pi<double>*light_speed<double>/
        reference_length;

    test_expected_pair_number <unit_system::SI, double>();
    test_expected_pair_number <unit_system::norm_omega, double>(reference_omega);
    test_expected_pair_number <unit_system::norm_lambda, double>(reference_length);
    test_expected_pair_number <unit_system::heaviside_lorentz, double>();
    test_expected_pair_number <unit_system::SI, float>();
    test_expected_pair_number <unit_system::norm_omega, float>(reference_omega);
    test_expected_pair_number <unit_system::norm_lambda, float>(reference_length);
    test_expected_pair_number <unit_system::heaviside_lorentz, float>();
}
