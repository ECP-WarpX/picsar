#ifndef PICSAR_MULTIPHYSICS_PHYS_CONSTANTS
#define PICSAR_MULTIPHYSICS_PHYS_CONSTANTS

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Uses some mathematical constants
#include "picsar_qed/math/math_constants.h"

namespace picsar::multi_physics::phys
{
    // Some useful physical constants in SI units (CODATA 2022)
    // (note that iIntermediate calculations, when required, are performed in
    // double precision to avoid numerical issues)
    template<typename RealType = double>
    constexpr auto electron_mass = RealType(9.1093837139e-31);

    template<typename RealType = double>
    constexpr auto elementary_charge = RealType(1.602176634e-19);

    template<typename RealType = double>
    constexpr auto light_speed = RealType(299792458.);

    template<typename RealType = double>
    constexpr auto reduced_plank = RealType(6.62607015e-34/(2.0*math::pi<>));

    template<typename RealType = double>
    constexpr auto vacuum_permittivity = RealType(8.8541878188e-12);

    // NOTE This is adjusted from the CODATA 2022 value 1.25663706127e-6,
    // so that the relation between exp_light_speed, exp_vacuum_permittivity,
    // and exp_vacuum_permeability is exact
    template<typename RealType = double>
    constexpr auto vacuum_permeability = RealType(1.2566370612685e-6);

    // NOTE This is calculated from alpha = mu_0/(4*pi)*q_e*q_e*c/hbar
    // and differs slightly from the CODATA 2022 value 0.0072973525643
    template<typename RealType = double>
    constexpr auto fine_structure =  RealType(0.0072973525643330135);

    template<typename RealType = double>
    constexpr auto eV = RealType(elementary_charge<>);

    template<typename RealType = double>
    constexpr auto KeV = RealType(elementary_charge<>*1e3);

    template<typename RealType = double>
    constexpr auto MeV = RealType(elementary_charge<>*1e6);

    template<typename RealType = double>
    constexpr auto GeV = RealType(elementary_charge<>*1e9);

    template<typename RealType = double>
    constexpr auto classical_electron_radius = RealType(2.8179403205e-15);

    //This constant is used for the Heaviside Lorentz unit system
    //(unfortunately, sqrt is not constexpr)
    template<typename RealType = double>
    constexpr auto sqrt_4_pi_fine_structure =
        RealType(0.3028221207690299);

    //
    template<typename RealType = double>
    constexpr auto schwinger_field = RealType(
        electron_mass<>*electron_mass<>*(light_speed<>*light_speed<>*light_speed<>)/
        (elementary_charge<>*reduced_plank<>));

    template<typename RealType = double>
    constexpr auto tau_e = RealType(classical_electron_radius<>/light_speed<>);
    //_______________________________________________________________________
}

#endif //PICSAR_MULTIPHYSICS_MATH_CONSTANTS
