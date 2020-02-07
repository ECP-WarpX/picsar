#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLE_GEN
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLE_GEN

//This .hpp file contais the implementation of the
//core functions of the nonlinear Breit Wheeler engine.
//If desired, these functions can be used directly
//without the higher-level interface.

#include <cmath>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar arrays
#include "../../math/vec_functions.hpp"

//Uses vector functions
#include "../../math/vec_functions.hpp"

//Uses chi functions
#include "../chi_functions.hpp"

//Uses physical constants
#include "../phys_constants.h"

//Uses unit conversion"
#include "../unit_conversion.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

 template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_aux_x(const RealType chi_phot, const RealType chi_ele)
    {
        return pow(
            chi_phot/(chi_ele*(chi_phot-chi_ele)),
            static_cast<RealType>(2.0/3.0));
    }

    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_aux_T_integrand(const RealType chi_phot, const RealType chi_ele)
    {
        if (chi_phot == static_cast<RealType>(0.0) ||
            chi_ele >= chi_phot ||
            chi_ele == static_cast<RealType>(0.0))
            return static_cast<RealType>(0.0);

        const auto xx = compute_x(chi_phot, chi_ele);
        const auto xx_3_2 = pow(xx, static_cast<RealType>(3.0/2.0));

        const auto div = static_cast<RealType>(math::pi*sqrt(3.0));

        const auto inner_integral = math::quad_a_inf<RealType>(
            [=](RealType s){
            /*  if (s >= static_cast<RealType>(__breit_wheeler_special_func_big_arg))
                return zero;*/
                return sqrt(s)*k_v(
                    static_cast<RealType>(1.0/3.0),
                    static_cast<RealType>(2.0/3.0)*pow(s, static_cast<RealType>(3.0/2.0)));
            }, xx);

        return (inner_integral-(static_cast<RealType>(2.0)-chi_phot*xx_3_2)
            *k_v(static_cast<RealType>(2.0/3.0),
                static_cast<RealType>(2.0/3.0)*xx_3_2))/div;
    }

    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_aux_T_function(const RealType chi_phot)
    {
        return math::quad_a_b<RealType>(
            compute_aux_T_integrand<RealType>,
            static_cast<RealType>(0.0), chi_phot);
    }

    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_pair_production_rate(
        const RealType t_energy_phot, const RealType chi_phot,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto energy_phot = t_energy_phot*
            fact_energy_to_SI_from<UnitSystem, RealType>();

        if(energy_phot ==  static_cast<RealType>(0.0) ||
            chi_phot ==  static_cast<RealType>(0.0))
            return static_cast<RealType>(0.0);

        const auto coeff = static_cast<RealType>(
            fine_structure * pow(electron_mass,2.0)*pow(light_speed,4.0)/
            reduced_plank)/( chi_phot * energy_phot);

        return coeff*compute_aux_T_function(chi_phot)*
            fact_rate_from_SI_to<RealType,UnitSystem>(lambda);
    }


}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLE_GEN
