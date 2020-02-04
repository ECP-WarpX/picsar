//####### Test module for unit conversion ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/unit_conversion"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "math_constants.h"
#include "phys_constants.h"
#include "unit_conversion.hpp"

using namespace picsar::multi_physics::phys;
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

//***Test mass conversion from and to SI

template<typename RealType>
void test_case_mass_from_SI()
{
    const auto mass = static_cast<RealType>(electron_mass);
    const auto res_SI = mass*fact_mass_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = mass*fact_mass_from_SI_to<unit_system::norm_omega, RealType>();
    const auto res_lambda = mass*fact_mass_from_SI_to<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(electron_mass);
    const auto exp_omega = static_cast<RealType>(1.0);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
void test_case_mass_to_SI()
{
    const auto mass = static_cast<RealType>(1.0);
    const auto res_SI = mass*fact_mass_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega = mass*fact_mass_to_SI_from<unit_system::norm_omega, RealType>();
    const auto res_lambda = mass*fact_mass_to_SI_from<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(electron_mass);
    const auto exp_lambda = static_cast<RealType>(electron_mass);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_mass )
{
    test_case_mass_from_SI<double>();
    test_case_mass_from_SI<float>();
    test_case_mass_to_SI<double>();
    test_case_mass_to_SI<float>();
}

//*******************************

//***Test charge conversion from and to SI

template<typename RealType>
constexpr void test_case_charge_from_SI()
{
    const auto charge = static_cast<RealType>(elementary_charge);
    const auto res_SI = charge*fact_charge_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = charge*fact_charge_from_SI_to<unit_system::norm_omega, RealType>();
    const auto res_lambda = charge*fact_charge_from_SI_to<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(elementary_charge);
    const auto exp_omega = static_cast<RealType>(1.0);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_charge_to_SI()
{
    const auto charge = static_cast<RealType>(1.0);;
    const auto res_SI = charge*fact_charge_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega = charge*fact_charge_to_SI_from<unit_system::norm_omega, RealType>();
    const auto res_lambda = charge*fact_charge_to_SI_from<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(elementary_charge);
    const auto exp_lambda = static_cast<RealType>(elementary_charge);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_charge )
{
    test_case_charge_from_SI<double>();
    test_case_charge_from_SI<float>();
    test_case_charge_to_SI<double>();
    test_case_charge_to_SI<float>();
}

//*******************************

//***Test velocity conversion from and to SI

template<typename RealType>
constexpr void test_case_velocity_from_SI()
{
    const auto velocity = static_cast<RealType>(light_speed);
    const auto res_SI = velocity*fact_velocity_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = velocity*fact_velocity_from_SI_to<unit_system::norm_omega, RealType>();
    const auto res_lambda = velocity*fact_velocity_from_SI_to<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(light_speed);
    const auto exp_omega = static_cast<RealType>(1.0);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_velocity_to_SI()
{
    const auto velocity = static_cast<RealType>(1.0);
    const auto res_SI = velocity*fact_velocity_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega = velocity*fact_velocity_to_SI_from<unit_system::norm_omega, RealType>();
    const auto res_lambda = velocity*fact_velocity_to_SI_from<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(light_speed);
    const auto exp_lambda = static_cast<RealType>(light_speed);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_velocity )
{
    test_case_velocity_from_SI<double>();
    test_case_velocity_from_SI<float>();
    test_case_velocity_to_SI<double>();
    test_case_velocity_to_SI<float>();
}

//*******************************

//***Test momentum conversion from and to SI

template<typename RealType>
constexpr void test_case_momentum_from_SI()
{
    const auto momentum = static_cast<RealType>(light_speed*electron_mass);
    const auto res_SI = momentum*fact_momentum_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = momentum*fact_momentum_from_SI_to<unit_system::norm_omega, RealType>();
    const auto res_lambda = momentum*fact_momentum_from_SI_to<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(light_speed*electron_mass);
    const auto exp_omega = static_cast<RealType>(1.0);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_momentum_to_SI()
{
    const auto momentum = static_cast<RealType>(1.0);
    const auto res_SI = momentum*fact_momentum_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega = momentum*fact_momentum_to_SI_from<unit_system::norm_omega, RealType>();
    const auto res_lambda = momentum*fact_momentum_to_SI_from<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(light_speed*electron_mass);
    const auto exp_lambda = static_cast<RealType>(light_speed*electron_mass);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_momentum )
{
    test_case_momentum_from_SI<double>();
    test_case_momentum_from_SI<float>();
    test_case_momentum_to_SI<double>();
    test_case_momentum_to_SI<float>();
}

//*******************************

//***Test length conversion from and to SI

template<typename RealType>
constexpr void test_case_length_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto length = static_cast<RealType>(lambda);
    const auto res_SI = length*
        fact_length_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = length*
        fact_length_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = length*
        fact_length_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(lambda);
    const auto exp_omega = static_cast<RealType>(2.0*pi);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_length_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto length = static_cast<RealType>(1.0);
    const auto res_SI =
        length*fact_length_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        length*fact_length_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        length*fact_length_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(lambda/2.0/pi);
    const auto exp_lambda = static_cast<RealType>(lambda);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_length )
{
    test_case_length_from_SI<double>();
    test_case_length_from_SI<float>();
    test_case_length_to_SI<double>();
    test_case_length_to_SI<float>();
}

//*******************************

//***Test area conversion from and to SI

template<typename RealType>
constexpr void test_case_area_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto area = static_cast<RealType>(lambda*lambda);
    const auto res_SI = area*
        fact_area_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = area*
        fact_area_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = area*
        fact_area_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(lambda*lambda);
    const auto exp_omega = static_cast<RealType>(2.0*pi*2.0*pi);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_area_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto area = static_cast<RealType>(1.0);
    const auto res_SI =
        area*fact_area_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        area*fact_area_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        area*fact_area_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(lambda*lambda/4.0/(pi*pi));
    const auto exp_lambda = static_cast<RealType>(lambda*lambda);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_area )
{
    test_case_area_from_SI<double>();
    test_case_area_from_SI<float>();
    test_case_area_to_SI<double>();
    test_case_area_to_SI<float>();
}

//*******************************

//***Test area conversion from and to SI

template<typename RealType>
constexpr void test_case_volume_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto volume = static_cast<RealType>(lambda*lambda*lambda);
    const auto res_SI = volume*
        fact_volume_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = volume*
        fact_volume_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = volume*
        fact_volume_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(lambda*lambda*lambda);
    const auto exp_omega = static_cast<RealType>(8.0*pi*pi*pi);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_volume_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto volume = static_cast<RealType>(1.0);
    const auto res_SI =
        volume*fact_volume_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        volume*fact_volume_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        volume*fact_volume_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(lambda*lambda*lambda/8.0/(pi*pi*pi));
    const auto exp_lambda = static_cast<RealType>(lambda*lambda*lambda);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_volume )
{
    test_case_volume_from_SI<double>();
    test_case_volume_from_SI<float>();
    test_case_volume_to_SI<double>();
    test_case_volume_to_SI<float>();
}

//*******************************

//***Test time conversion from and to SI

template<typename RealType>
constexpr void test_case_time_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto time = static_cast<RealType>(lambda/light_speed);
    const auto res_SI = time*
        fact_time_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = time*
        fact_time_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = time*
        fact_time_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(lambda/light_speed);
    const auto exp_omega = static_cast<RealType>(2.0*pi);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_time_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto time = static_cast<RealType>(1.0);
    const auto res_SI =
        time*fact_time_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        time*fact_time_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        time*fact_time_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(lambda/2.0/pi/light_speed);
    const auto exp_lambda = static_cast<RealType>(lambda/light_speed);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_time )
{
    test_case_time_from_SI<double>();
    test_case_time_from_SI<float>();
    test_case_time_to_SI<double>();
    test_case_time_to_SI<float>();
}

//*******************************

//***Test rate conversion from and to SI

template<typename RealType>
constexpr void test_case_rate_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto rate = static_cast<RealType>(light_speed/lambda);
    const auto res_SI = rate*
        fact_rate_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = rate*
        fact_rate_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = rate*
        fact_rate_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(light_speed/lambda);
    const auto exp_omega = static_cast<RealType>(1.0/2.0/pi);
    const auto exp_lambda = static_cast<RealType>(1.0);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_rate_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto rate = static_cast<RealType>(1.0);
    const auto res_SI =
        rate*fact_rate_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        rate*fact_rate_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        rate*fact_rate_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(2.0*pi*light_speed/lambda);
    const auto exp_lambda = static_cast<RealType>(light_speed/lambda);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_rate )
{
    test_case_rate_from_SI<double>();
    test_case_rate_from_SI<float>();
    test_case_rate_to_SI<double>();
    test_case_rate_to_SI<float>();
}

//*******************************

//***Test E conversion from and to SI

template<typename RealType>
constexpr void test_case_E_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto E = static_cast<RealType>(
        electron_mass*light_speed*2.0*pi*light_speed/elementary_charge)/lambda;
    const auto res_SI = E*
        fact_E_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = E*
        fact_E_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = E*
        fact_E_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(E);
    const auto exp_omega = static_cast<RealType>(1.0);
    const auto exp_lambda = static_cast<RealType>(2.0*pi);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_E_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto E = static_cast<RealType>(1.0);
    const auto res_SI =
        E*fact_E_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        E*fact_E_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        E*fact_E_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(
        electron_mass*light_speed*2.0*pi*light_speed/elementary_charge)/lambda;
    const auto exp_lambda = static_cast<RealType>(
        electron_mass*light_speed*light_speed/elementary_charge)/lambda;

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_E )
{
    test_case_E_from_SI<double>();
    test_case_E_from_SI<float>();
    test_case_E_to_SI<double>();
    test_case_E_to_SI<float>();
}

//*******************************

//***Test B conversion from and to SI

template<typename RealType>
constexpr void test_case_B_from_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto B = static_cast<RealType>(
        electron_mass*light_speed*2.0*pi/elementary_charge)/lambda;
    const auto res_SI = B*
        fact_B_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = B*
        fact_B_from_SI_to<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda = B*
        fact_B_from_SI_to<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(B);
    const auto exp_omega = static_cast<RealType>(1.0);
    const auto exp_lambda = static_cast<RealType>(2.0*pi);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_B_to_SI()
{
    const auto lambda = static_cast<RealType>(800.0e-9);
    const auto omega = static_cast<RealType>(2.0*pi*light_speed/lambda);

    const auto B = static_cast<RealType>(1.0);
    const auto res_SI =
        B*fact_B_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        B*fact_B_to_SI_from<unit_system::norm_omega, RealType>(omega);
    const auto res_lambda =
        B*fact_B_to_SI_from<unit_system::norm_lambda, RealType>(lambda);
    const auto exp_SI = static_cast<RealType>(1.0);
    const auto exp_omega = static_cast<RealType>(
        electron_mass*light_speed*2.0*pi/elementary_charge)/lambda;
    const auto exp_lambda = static_cast<RealType>(
        electron_mass*light_speed/elementary_charge)/lambda;

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_B )
{
    test_case_B_from_SI<double>();
    test_case_B_from_SI<float>();
    test_case_B_to_SI<double>();
    test_case_B_to_SI<float>();
}

//*******************************

//***Test energy conversion from and to SI

template<typename RealType>
constexpr void test_case_energy_from_SI()
{
    const auto energy = static_cast<RealType>(
        electron_mass*light_speed*light_speed);
    const auto res_SI = energy*
        fact_energy_from_SI_to<unit_system::SI, RealType>();
    const auto res_omega = energy*
        fact_energy_from_SI_to<unit_system::norm_omega, RealType>();
    const auto res_lambda = energy*
        fact_energy_from_SI_to<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(energy);
    const auto exp_omega = static_cast<RealType>(energy/
        (electron_mass*light_speed*light_speed));
    const auto exp_lambda = static_cast<RealType>(energy/
        (electron_mass*light_speed*light_speed));

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

template<typename RealType>
constexpr void test_case_energy_to_SI()
{
    const auto energy = static_cast<RealType>(1.0);
    const auto res_SI =
        energy*fact_energy_to_SI_from<unit_system::SI, RealType>();
    const auto res_omega =
        energy*fact_energy_to_SI_from<unit_system::norm_omega, RealType>();
    const auto res_lambda =
        energy*fact_energy_to_SI_from<unit_system::norm_lambda, RealType>();
    const auto exp_SI = static_cast<RealType>(energy);
    const auto exp_omega = static_cast<RealType>(
        energy*electron_mass*light_speed*light_speed);
    const auto exp_lambda = static_cast<RealType>(
        energy*electron_mass*light_speed*light_speed);

    BOOST_CHECK_SMALL((res_SI-exp_SI)/exp_SI, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_omega-exp_omega)/exp_omega, tolerance<RealType>());
    BOOST_CHECK_SMALL((res_lambda-exp_lambda)/exp_lambda, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_energy )
{
    test_case_energy_from_SI<double>();
    test_case_energy_from_SI<float>();
    test_case_energy_to_SI<double>();
    test_case_energy_to_SI<float>();
}

//*******************************
