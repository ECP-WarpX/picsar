#ifndef PICSAR_MULTIPHYSICS_PHYS_CONSTANTS
#define PICSAR_MULTIPHYSICS_PHYS_CONSTANTS

#include "../qed_commons.h"
#include "../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace phys{

    // Physical constants in SI units

    template<typename RealType = double>
    constexpr auto electron_mass = RealType(9.10938356e-31);

    template<typename RealType = double>
    constexpr auto elementary_charge = RealType(1.6021766208e-19);

    template<typename RealType = double>
    constexpr auto light_speed = RealType(299792458.);

    template<typename RealType = double>
    constexpr auto reduced_plank = RealType(1.054571800e-34);

    template<typename RealType = double>
    constexpr auto vacuum_permittivity = RealType(8.854187817e-12);

    template<typename RealType = double>
    constexpr auto fine_structure =  RealType(0.0072973525664);

    //Intermediate calculations of the following quantities are performed with
    //double precision to avoid numerical issues
    template<typename RealType = double>
    constexpr auto classical_electron_radius = RealType(
        elementary_charge*elementary_charge /
        (4.0*math::pi * vacuum_permittivity *
        electron_mass * light_speed * light_speed));

    template<typename RealType = double>
    constexpr auto schwinger_field = RealType(
        electron_mass*electron_mass*(light_speed*light_speed*light_speed)/
        (elementary_charge*reduced_plank));

    template<typename RealType = double>
    constexpr auto tau_e = RealType(classical_electron_radius/light_speed);
    //_______________________________________________________________________
}
}
}

#endif //PICSAR_MULTIPHYSICS_MATH_CONSTANTS
