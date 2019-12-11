#ifndef PICSAR_MULTIPHYSICS_PHYS_CONSTANTS
#define PICSAR_MULTIPHYSICS_PHYS_CONSTANTS

#include "../qed_commons.h"
#include "../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace phys{
    constexpr const double electron_mass = 9.10938356e-31;
    constexpr const double elementary_charge = 1.6021766208e-19;
    constexpr const double light_speed = 299792458.;
    constexpr const double reduced_plank = 1.054571800e-34;
    constexpr const double vacuum_permittivity =  8.854187817e-12;
    constexpr const double fine_structure =  0.0072973525664;

    constexpr const double classical_electron_radius =
        elementary_charge*elementary_charge /
        (4.0*math::pi*vacuum_permittivity*electron_mass*light_speed*light_speed);

    constexpr const double schwinger_field =
        electron_mass*electron_mass*(light_speed*light_speed*light_speed)/
        (elementary_charge*reduced_plank);

    constexpr const double pair_prod_rate_coeff =
        fine_structure * electron_mass * light_speed * light_speed /
        (reduced_plank);

    constexpr const double tau_e = classical_electron_radius/light_speed;

    constexpr const double quantum_synchrotron_rate_coeff =
        fine_structure*fine_structure/tau_e;

    constexpr const double schwinger_pair_prod_coeff =
          elementary_charge*elementary_charge*
          schwinger_field*schwinger_field/
          (4.0*math::pi*math::pi*reduced_plank*reduced_plank*light_speed);
}
}
}

#endif //PICSAR_MULTIPHYSICS_MATH_CONSTANTS