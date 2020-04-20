//####### Test module for unit conversion ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/unit_conversion"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

#include <array>
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

enum quantity_no_refs{
    mass,
    charge,
    velocity,
    momentum,
    energy
};

enum quantity_w_refs{
    length,
    area,
    volume,
    ttime,
    rate,
    E,
    B
};

template <quantity_no_refs Qnt>
struct fact_no_ref_to_SI{
    template<unit_system From, typename RealType>
    static constexpr RealType from();
};

template<>
struct fact_no_ref_to_SI<quantity_no_refs::mass>
{
    template<unit_system From, typename RealType>
    static constexpr RealType from(){
        return fact_mass_to_SI_from<From, RealType>();}
};

template<>
struct fact_no_ref_to_SI<quantity_no_refs::charge>
{
    template<unit_system From, typename RealType>
    static constexpr RealType from(){
        return fact_charge_to_SI_from<From, RealType>();}
};

template<>
struct fact_no_ref_to_SI<quantity_no_refs::velocity>
{
    template<unit_system From, typename RealType>
    static constexpr RealType from(){
        return fact_velocity_to_SI_from<From, RealType>();}
};

template<>
struct fact_no_ref_to_SI<quantity_no_refs::momentum>
{
    template<unit_system From, typename RealType>
    static constexpr RealType from(){
        return fact_momentum_to_SI_from<From, RealType>();}
};

template<>
struct fact_no_ref_to_SI<quantity_no_refs::energy>
{
    template<unit_system From, typename RealType>
    static constexpr RealType from(){
        return fact_energy_to_SI_from<From, RealType>();}
};


template<typename RealType, quantity_no_refs Quantity>
constexpr void test_to_SI_no_refs(val_pack<RealType> vals)
{
    const auto fact_SI = fact_no_ref_to_SI<Quantity>::
        template from<unit_system::SI, RealType>();
    const auto fact_omega = fact_no_ref_to_SI<Quantity>::
        template from<unit_system::norm_omega, RealType>();
    const auto fact_lambda = fact_no_ref_to_SI<Quantity>::
        template from<unit_system::norm_lambda, RealType>();
    const auto fact_hl = fact_no_ref_to_SI<Quantity>::
        template from<unit_system::heaviside_lorentz, RealType>();

    const auto res_SI2SI = vals.SI*fact_SI;
    const auto res_omega2SI = vals.omega*fact_omega;
    const auto res_lambda2SI = vals.lambda*fact_lambda;
    const auto res_hl2SI = vals.hl*fact_hl;
    const auto all_res = std::array<RealType,4>{
        res_SI2SI, res_omega2SI, res_lambda2SI, res_hl2SI};
    for (const auto& res : all_res)
        BOOST_CHECK_SMALL((res-vals.SI)/vals.SI, tolerance<RealType>());
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

// ***Test mass conversion to SI
template<typename RealType>
void test_case_mass_to_SI()
{
    constexpr auto mass_SI = electron_mass<RealType>;
    constexpr auto mass_omega = static_cast<RealType>(1.0);
    constexpr auto mass_lambda = static_cast<RealType>(1.0);
    constexpr auto mass_hl = electron_mass<RealType>*light_speed<RealType>*
        light_speed<RealType>/heaviside_lorentz_reference_energy<RealType>;
    constexpr auto all_masses = val_pack<RealType>{mass_SI, mass_omega, mass_lambda, mass_hl};

    test_to_SI_no_refs<RealType, quantity_no_refs::mass>(all_masses);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_mass_to_SI )
{
    test_case_mass_to_SI<double>();
    test_case_mass_to_SI<float>();
}

// ***Test charge conversion to SI
template<typename RealType>
void test_case_charge_to_SI()
{
    constexpr auto charge_SI = elementary_charge<RealType>;
    constexpr auto charge_omega = static_cast<RealType>(1.0);
    constexpr auto charge_lambda = static_cast<RealType>(1.0);
    constexpr auto charge_hl = sqrt_4_pi_fine_structure<RealType>;
    constexpr auto all_charges = val_pack<RealType>{charge_SI, charge_omega, charge_lambda, charge_hl};

    test_to_SI_no_refs<RealType, quantity_no_refs::charge>(all_charges);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_charge_to_SI )
{
    test_case_charge_to_SI<double>();
    test_case_charge_to_SI<float>();
}

// ***Test velocity conversion to SI
template<typename RealType>
void test_case_velocity_to_SI()
{
    constexpr auto velocity_SI = light_speed<RealType>;
    constexpr auto velocity_omega = static_cast<RealType>(1.0);
    constexpr auto velocity_lambda = static_cast<RealType>(1.0);
    constexpr auto velocity_hl = static_cast<RealType>(1.0);
    constexpr auto all_velocities = val_pack<RealType>{velocity_SI, velocity_omega, velocity_lambda, velocity_hl};

    test_to_SI_no_refs<RealType, quantity_no_refs::velocity>(all_velocities);
}

BOOST_AUTO_TEST_CASE( picsar_unit_conv_velocity_to_SI )
{
    test_case_velocity_to_SI<double>();
    test_case_velocity_to_SI<float>();
}
