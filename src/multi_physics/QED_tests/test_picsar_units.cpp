//####### Test module for unit conversion ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/unit_conversion"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

#include <array>
#include <algorithm>
#include <functional>

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

//Auxiliary functions for tests
template<typename RealType>
struct val_pack{
    RealType SI;
    RealType omega;
    RealType lambda;
    RealType hl;
};

template<typename RealType, quantity Quantity>
constexpr void test_to_SI(val_pack<RealType> vals,
    RealType reference_omega = -1.0,
    RealType reference_length = -1.0)
{
    const auto fact_SI =
        conv<Quantity, unit_system::SI,
            unit_system::SI, RealType>::fact();
    const auto fact_omega =
        conv<Quantity, unit_system::norm_omega,
            unit_system::SI, RealType>::fact(reference_omega);
    const auto fact_lambda =
        conv<Quantity, unit_system::norm_lambda,
            unit_system::SI, RealType>::fact(reference_length);
    const auto fact_hl =
        conv<Quantity, unit_system::heaviside_lorentz,
            unit_system::SI, RealType>::fact();

    const auto res_SI2SI = vals.SI*fact_SI;
    const auto res_omega2SI = vals.omega*fact_omega;
    const auto res_lambda2SI = vals.lambda*fact_lambda;
    const auto res_hl2SI = vals.hl*fact_hl;
    const auto all_res = std::array<RealType,4>{
        res_SI2SI, res_omega2SI, res_lambda2SI, res_hl2SI};
    for (const auto& res : all_res)
        BOOST_CHECK_SMALL((res-vals.SI)/vals.SI, tolerance<RealType>());
}

template<typename RealType, quantity Quantity>
constexpr void test_from_to(
    val_pack<RealType> vals,
    RealType reference_omega = -1.0,
    RealType reference_length = -1.0)

{
    const auto from_SI_to_all = std::array<RealType,4>
    {
        vals.SI*conv<Quantity, unit_system::SI,
            unit_system::SI, RealType>::fact(),
        vals.SI*conv<Quantity, unit_system::SI,
            unit_system::norm_omega, RealType>::fact(1.0,reference_omega),
        vals.SI*conv<Quantity, unit_system::SI,
            unit_system::norm_lambda, RealType>::fact(1.0,reference_length),
        vals.SI*conv<Quantity, unit_system::SI,
            unit_system::heaviside_lorentz, RealType>::fact(),
    };

    const auto from_omega_to_all = std::array<RealType,4>
    {
        vals.omega*conv<Quantity, unit_system::norm_omega,
            unit_system::SI, RealType>::fact(reference_omega),
        vals.omega*conv<Quantity, unit_system::norm_omega,
            unit_system::norm_omega, RealType>::fact(reference_omega, reference_omega),
        vals.omega*conv<Quantity, unit_system::norm_omega,
            unit_system::norm_lambda, RealType>::fact(reference_omega, reference_length),
        vals.omega*conv<Quantity, unit_system::norm_omega,
            unit_system::heaviside_lorentz, RealType>::fact(reference_omega),
    };

    const auto from_lambda_to_all = std::array<RealType,4>
    {
        vals.lambda*conv<Quantity, unit_system::norm_lambda,
            unit_system::SI, RealType>::fact(reference_length),
        vals.lambda*conv<Quantity, unit_system::norm_lambda,
            unit_system::norm_omega, RealType>::fact(reference_length, reference_omega),
        vals.lambda*conv<Quantity, unit_system::norm_lambda,
            unit_system::norm_lambda, RealType>::fact(reference_length, reference_length),
        vals.lambda*conv<Quantity, unit_system::norm_lambda,
            unit_system::heaviside_lorentz, RealType>::fact(reference_length, 1.0),
    };

    const auto from_hl_to_all = std::array<RealType,4>
    {
        vals.hl*conv<Quantity, unit_system::heaviside_lorentz,
            unit_system::SI, RealType>::fact(),
        vals.hl*conv<Quantity, unit_system::heaviside_lorentz,
            unit_system::norm_omega, RealType>::fact(1.0,reference_omega),
        vals.hl*conv<Quantity, unit_system::heaviside_lorentz,
            unit_system::norm_lambda, RealType>::fact(1.0,reference_length),
        vals.hl*conv<Quantity, unit_system::heaviside_lorentz,
            unit_system::heaviside_lorentz, RealType>::fact(),
    };

    const auto fact_SI = conv<Quantity, unit_system::SI,
        unit_system::SI, RealType>::fact();
    const auto fact_omega = conv<Quantity, unit_system::norm_omega,
        unit_system::SI, RealType>::fact(reference_omega);
    const auto fact_lambda = conv<Quantity, unit_system::norm_lambda,
        unit_system::SI, RealType>::fact(reference_length);
    const auto fact_hl = conv<Quantity, unit_system::heaviside_lorentz,
        unit_system::SI, RealType>::fact();

    const auto all_facts = std::array<RealType, 4>{
        fact_SI, fact_omega, fact_lambda, fact_hl};

    const auto all_data = std::array<std::array<RealType, 4>, 4>{
        from_SI_to_all,
        from_omega_to_all,
        from_lambda_to_all,
        from_hl_to_all
    };

    for (auto data: all_data)
    {
        std::transform( data.begin(), data.end(),
            all_facts.begin(), data.begin(),
            std::multiplies<RealType>());

        for (const auto& res : data){
            BOOST_CHECK_SMALL((res-vals.SI)/vals.SI, tolerance<RealType>());
        }
    }
}

// ------------- Tests --------------

// ***Test energy reference for Heaviside Lorentz units
template<typename RealType>
void test_case_hl_reference_energy()
{
    BOOST_CHECK_SMALL(
        (heaviside_lorentz_reference_energy<RealType>-MeV<RealType>)/MeV<RealType>,
        tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_heaviside_lorentz_ref_energy )
{
    test_case_hl_reference_energy<double>();
    test_case_hl_reference_energy<float>();
}

// ***Test electron rest energy in Heaviside Lorentz units
template<typename RealType>
void test_case_hl_electron_rest_energy()
{
    constexpr auto exp = static_cast<RealType>(
        electron_mass<double>*light_speed<double>*light_speed<double>/
        MeV<double>);
    constexpr auto res = heaviside_lorentz_electron_rest_energy<RealType>;
    BOOST_CHECK_SMALL((res-exp)/exp, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_heaviside_lorentz_electron_rest_energy )
{
    test_case_hl_electron_rest_energy<double>();
    test_case_hl_electron_rest_energy<float>();
}

// ***Test Schwinger field in Heaviside Lorentz units
template<typename RealType>
void test_case_hl_schwinger_field()
{
    constexpr auto exp = static_cast<RealType>(
        schwinger_field<double>*conv<quantity::E,
            unit_system::SI, unit_system::heaviside_lorentz, double>::fact());
    constexpr auto res = heaviside_lorentz_schwinger_field<RealType>;
    BOOST_CHECK_SMALL((res-exp)/exp, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_heaviside_lorentz_schwinger_field )
{
    test_case_hl_schwinger_field<double>();
    test_case_hl_schwinger_field<float>();
}

// ***Test mass conversion to SI and all to all
template<typename RealType>
void test_case_mass()
{
    constexpr auto mass_SI = electron_mass<RealType>;
    constexpr auto mass_omega = static_cast<RealType>(1.0);
    constexpr auto mass_lambda = static_cast<RealType>(1.0);
    constexpr auto mass_hl = electron_mass<RealType>*light_speed<RealType>*
        light_speed<RealType>/heaviside_lorentz_reference_energy<RealType>;
    constexpr auto all_masses = val_pack<RealType>{mass_SI, mass_omega, mass_lambda, mass_hl};

    test_to_SI<RealType, quantity::mass>(all_masses);
    test_from_to<RealType, quantity::mass>(all_masses);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_mass )
{
    test_case_mass<double>();
    test_case_mass<float>();
}

// ***Test charge conversion to SI and all to all
template<typename RealType>
void test_case_charge()
{
    constexpr auto charge_SI = elementary_charge<RealType>;
    constexpr auto charge_omega = static_cast<RealType>(1.0);
    constexpr auto charge_lambda = static_cast<RealType>(1.0);
    constexpr auto charge_hl = sqrt_4_pi_fine_structure<RealType>;
    constexpr auto all_charges = val_pack<RealType>{charge_SI, charge_omega, charge_lambda, charge_hl};

    test_to_SI<RealType, quantity::charge>(all_charges);
    test_from_to<RealType, quantity::charge>(all_charges);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_charge )
{
    test_case_charge<double>();
    test_case_charge<float>();
}

// ***Test velocity conversion to SI and all to all
template<typename RealType>
void test_case_velocity()
{
    constexpr auto velocity_SI = light_speed<RealType>;
    constexpr auto velocity_omega = static_cast<RealType>(1.0);
    constexpr auto velocity_lambda = static_cast<RealType>(1.0);
    constexpr auto velocity_hl = static_cast<RealType>(1.0);
    constexpr auto all_velocities = val_pack<RealType>{velocity_SI, velocity_omega, velocity_lambda, velocity_hl};

    test_to_SI<RealType, quantity::velocity>(all_velocities);
    test_from_to<RealType, quantity::velocity>(all_velocities);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_velocity )
{
    test_case_velocity<double>();
    test_case_velocity<float>();
}

// ***Test momentum conversion to SI and all to all
template<typename RealType>
void test_case_momentum()
{
    constexpr auto momentum_SI = electron_mass<RealType>*light_speed<RealType>;
    constexpr auto momentum_omega = static_cast<RealType>(1.0);
    constexpr auto momentum_lambda = static_cast<RealType>(1.0);
    constexpr auto momentum_hl = static_cast<RealType>(
        electron_mass<double>*light_speed<double>*light_speed<double>/
        heaviside_lorentz_reference_energy<double>);
    constexpr auto all_momenta = val_pack<RealType>{momentum_SI, momentum_omega, momentum_lambda, momentum_hl};

    test_to_SI<RealType, quantity::momentum>(all_momenta);
    test_from_to<RealType, quantity::momentum>(all_momenta);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_momentum )
{
    test_case_momentum<double>();
    test_case_momentum<float>();
}

// ***Test energy conversion to SI and all to all
template<typename RealType>
void test_case_energy()
{
    constexpr auto energy_SI = GeV<RealType>;
    constexpr auto energy_omega = static_cast<RealType>(
        GeV<double>/electron_mass<double>/light_speed<double>/light_speed<double>);
    constexpr auto energy_lambda = static_cast<RealType>(
        GeV<double>/electron_mass<double>/light_speed<double>/light_speed<double>);
    constexpr auto energy_hl = static_cast<RealType>(
        GeV<double>/heaviside_lorentz_reference_energy<double>);
    constexpr auto all_energies = val_pack<RealType>{energy_SI, energy_omega, energy_lambda, energy_hl};

    test_to_SI<RealType, quantity::energy>(all_energies);
    test_from_to<RealType, quantity::energy>(all_energies);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_energy )
{
    test_case_energy<double>();
    test_case_energy<float>();
}


// ***Test length conversion to SI and all to all
template<typename RealType>
void test_case_length()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr auto length_SI = reference_length;
    constexpr auto length_omega = static_cast<RealType>(2.0* pi<double>);
    constexpr auto length_lambda = static_cast<RealType>(1.0);
    constexpr auto length_hl = static_cast<RealType>(
        heaviside_lorentz_reference_energy<double>*reference_length/
        reduced_plank<double>/light_speed<double>);

    constexpr auto all_lenghts = val_pack<RealType>{length_SI, length_omega, length_lambda, length_hl};

    test_to_SI<RealType, quantity::length>(
        all_lenghts, reference_omega, reference_length);
    test_from_to<RealType, quantity::length>(
        all_lenghts, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_length)
{
    test_case_length<double>();
    test_case_length<float>();
}

// ***Test area conversion to SI and all to all
template<typename RealType>
void test_case_area_to_SI()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr auto area_SI = reference_length*reference_length;
    constexpr auto area_omega = static_cast<RealType>(4.0* pi<double>*pi<double>);
    constexpr auto area_lambda = static_cast<RealType>(1.0);
    constexpr auto area_hl = static_cast<RealType>(
        heaviside_lorentz_reference_energy<double>*
        heaviside_lorentz_reference_energy<double>*
        reference_length*reference_length/
        reduced_plank<double>/reduced_plank<double>/
        light_speed<double>/light_speed<double>);

    constexpr auto all_areas = val_pack<RealType>{area_SI, area_omega, area_lambda, area_hl};

    test_to_SI<RealType, quantity::area>(
        all_areas, reference_omega, reference_length);
    test_from_to<RealType, quantity::area>(
        all_areas, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_area_to_SI )
{
    test_case_area_to_SI<double>();
    test_case_area_to_SI<float>();
}

// ***Test volume conversion to SI and all to all
template<typename RealType>
void test_case_volume_to_SI()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr auto volume_SI = reference_length*reference_length*reference_length;
    constexpr auto volume_omega = static_cast<RealType>(8.0*
        pi<double>*pi<double>*pi<double>);
    constexpr auto volume_lambda = static_cast<RealType>(1.0);
    constexpr auto volume_hl = static_cast<RealType>(
        heaviside_lorentz_reference_energy<double>*
        heaviside_lorentz_reference_energy<double>*
        heaviside_lorentz_reference_energy<double>*
        reference_length*reference_length*reference_length/
        reduced_plank<double>/reduced_plank<double>/reduced_plank<double>/
        light_speed<double>/light_speed<double>/light_speed<double>);

    constexpr auto all_volumes = val_pack<RealType>{volume_SI, volume_omega, volume_lambda, volume_hl};

    test_to_SI<RealType, quantity::volume>(
        all_volumes, reference_omega, reference_length);
    test_from_to<RealType, quantity::volume>(
        all_volumes, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_volume_to_SI )
{
    test_case_volume_to_SI<double>();
    test_case_volume_to_SI<float>();
}

// ***Test time conversion to SI and all to all
template<typename RealType>
void test_case_time_to_SI()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr auto time_SI = static_cast<RealType>(reference_length/light_speed<double>);
    constexpr auto time_omega = static_cast<RealType>(2.0*pi<double>);
    constexpr auto time_lambda = static_cast<RealType>(1.0);
    constexpr auto time_hl = static_cast<RealType>(
        (reference_length/light_speed<double>)*
        heaviside_lorentz_reference_energy<double>/
        reduced_plank<double>);

    constexpr auto all_times = val_pack<RealType>{time_SI, time_omega, time_lambda, time_hl};

    test_to_SI<RealType, quantity::time>(
        all_times, reference_omega, reference_length);
    test_from_to<RealType, quantity::time>(
        all_times, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_time_to_SI )
{
    test_case_time_to_SI<double>();
    test_case_time_to_SI<float>();
}

// ***Test rate conversion to SI and all to all
template<typename RealType>
void test_case_rate_to_SI()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr auto rate_SI = static_cast<RealType>(light_speed<double>/reference_length);
    constexpr auto rate_omega = static_cast<RealType>(1/(2.0*pi<double>));
    constexpr auto rate_lambda = static_cast<RealType>(1.0);
    constexpr auto rate_hl = static_cast<RealType>(
        (light_speed<double>/reference_length)*
        reduced_plank<double>/
        heaviside_lorentz_reference_energy<double>);

    constexpr auto all_rates = val_pack<RealType>{rate_SI, rate_omega, rate_lambda, rate_hl};

    test_to_SI<RealType, quantity::rate>(
        all_rates, reference_omega, reference_length);
    test_from_to<RealType, quantity::rate>(
        all_rates, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_rate_to_SI )
{
    test_case_rate_to_SI<double>();
    test_case_rate_to_SI<float>();
}

// ***Test E conversion to SI and all to all
template<typename RealType>
void test_case_E_to_SI()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr double ref_field_omega = reference_omega*
        (electron_mass<double>*light_speed<double>/elementary_charge<double>);
    constexpr double ref_field_length = electron_mass<double>*
        light_speed<double>*light_speed<double>/elementary_charge<double>/reference_length;
    constexpr double ref_field_hl =
        heaviside_lorentz_reference_energy<double>*
        heaviside_lorentz_reference_energy<double>*
        sqrt_4_pi_fine_structure<double>/
        elementary_charge<double>/reduced_plank<double>/light_speed<double>;

    constexpr auto E_SI = schwinger_field<RealType>;
    constexpr auto E_omega = static_cast<RealType>(
        schwinger_field<double>/ref_field_omega);
    constexpr auto E_lambda = static_cast<RealType>(schwinger_field<double>/
        ref_field_length);
    constexpr auto E_hl = static_cast<RealType>(schwinger_field<double>/
        ref_field_hl);

    constexpr auto all_E = val_pack<RealType>{E_SI, E_omega, E_lambda, E_hl};

    test_to_SI<RealType, quantity::E>(
        all_E, reference_omega, reference_length);
    test_from_to<RealType, quantity::E>(
        all_E, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_E_to_SI )
{
    test_case_E_to_SI<double>();
    test_case_E_to_SI<float>();
}

// ***Test B conversion to SI and all to all
template<typename RealType>
void test_case_B_to_SI()
{
    constexpr auto reference_length = static_cast<RealType>(800.0e-9);
    constexpr auto reference_omega = static_cast<RealType>(
        2.0*pi<double>*light_speed<double>/reference_length);

    constexpr double ref_field_omega = reference_omega*
        (electron_mass<double>/elementary_charge<double>);
    constexpr double ref_field_length = electron_mass<double>*
        light_speed<double>/elementary_charge<double>/reference_length;
    constexpr double ref_field_hl =
        heaviside_lorentz_reference_energy<double>*
        heaviside_lorentz_reference_energy<double>*
        sqrt_4_pi_fine_structure<double>/
        elementary_charge<double>/reduced_plank<double>/
        light_speed<double>/light_speed<double>;

    constexpr double mag_schwinger = schwinger_field<double>/light_speed<double>;

    constexpr auto B_SI = static_cast<RealType>(mag_schwinger);
    constexpr auto B_omega = static_cast<RealType>(
        mag_schwinger/ref_field_omega);
    constexpr auto B_lambda = static_cast<RealType>(mag_schwinger/
        ref_field_length);
    constexpr auto B_hl = static_cast<RealType>(mag_schwinger/
        ref_field_hl);

    constexpr auto all_B = val_pack<RealType>{B_SI, B_omega, B_lambda, B_hl};

    test_to_SI<RealType, quantity::B>(
        all_B, reference_omega, reference_length);
    test_from_to<RealType, quantity::B>(
        all_B, reference_omega, reference_length);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_B_to_SI )
{
    test_case_B_to_SI<double>();
    test_case_B_to_SI<float>();
}
