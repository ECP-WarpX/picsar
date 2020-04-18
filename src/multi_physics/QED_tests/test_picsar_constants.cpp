//####### Test module for physical and mathematical constants ##################

//Define Module name
 #define BOOST_TEST_MODULE "phys/unit_conversion"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "math_constants.h"
#include "phys_constants.h"

using namespace picsar::multi_physics::phys;
using namespace picsar::multi_physics::math;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-14;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-7;

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

//***Test math constants

template<typename RealType>
void test_case_const_math()
{
    const auto exp_pi =
        static_cast<RealType>(3.14159265358979323846264338327950288);

    BOOST_CHECK_SMALL((pi<RealType>-exp_pi)/exp_pi, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_const_math )
{
    test_case_const_math<double>();
    test_case_const_math<float>();
}

//*******************************

//***Test physical constants

template<typename RealType>
void test_case_const_phys()
{
    const auto exp_electron_mass =
        static_cast<RealType>(9.1093837015e-31);
    const auto exp_elementary_charge =
        static_cast<RealType>(1.602176634e-19);
    const auto exp_light_speed =
        static_cast<RealType>(299792458.);
    const auto exp_reduced_plank =
        static_cast<RealType>(1.054571817e-34);
    const auto exp_vacuum_permittivity =
        static_cast<RealType>(8.8541878128e-12);
    const auto exp_fine_structure =
        static_cast<RealType>(0.0072973525693);
    const auto exp_classical_electron_radius =
        static_cast<RealType>(2.81794032620493e-15);
    const auto exp_schwinger_field =
        static_cast<RealType>(1.32328547494817e18);
    const auto exp_tau_e =
        static_cast<RealType>(9.39963715232933e-24);

    BOOST_CHECK_SMALL(
        (electron_mass<RealType>-exp_electron_mass)/exp_electron_mass,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (elementary_charge<RealType>-exp_elementary_charge)/exp_elementary_charge,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (light_speed<RealType>-exp_light_speed)/exp_light_speed,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (reduced_plank<RealType>-exp_reduced_plank)/exp_reduced_plank,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (vacuum_permittivity<RealType>-exp_vacuum_permittivity)/exp_vacuum_permittivity,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (fine_structure<RealType>-exp_fine_structure)/exp_fine_structure,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (classical_electron_radius<RealType>-exp_classical_electron_radius)/exp_classical_electron_radius,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (schwinger_field<RealType>-exp_schwinger_field)/exp_schwinger_field,
        tolerance<RealType>());

    BOOST_CHECK_SMALL(
        (tau_e<RealType>-exp_tau_e)/exp_tau_e,
        tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_const_phys )
{
    test_case_const_phys<double>();
    test_case_const_phys<float>();
}

//*******************************
