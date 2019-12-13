#ifndef PICSAR_MULTIPHYSICS_CHI_FUNCTIONS
#define PICSAR_MULTIPHYSICS_CHI_FUNCTIONS

//This .hpp file contains the functions to calcuate the chi parameter for
//photons and leptons

#include <cmath>

//Should be included by all the src files of the library
#include "../qed_commons.h"

//Uses picsar_arrays
#include "../containers/picsar_array.hpp"

#include "../math/vec_functions.hpp"

#include "unit_conversion.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{

    //chi for photons
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_photon(
        const math::vec3<RealType> t_p,
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType lambda = static_cast<RealType>(1.0))
    {

        const auto p = t_p * fact_momentum_to_SI_from<UnitSystem, RealType>();
        const auto em_e = t_em_e * fact_E_to_SI_from<UnitSystem, RealType>(lambda);
        const auto em_b = t_em_b * fact_B_to_SI_from<UnitSystem, RealType>(lambda);

        auto norm_p = math::norm(p);
        if(norm_p == static_cast<RealType>(0.0))
            return static_cast<RealType>(0.0);

        const auto p_unit = p / norm_p;
        const auto em_eperp = em_e - math::dot(p_unit,em_e)*p_unit;
        const auto mod = math::norm(em_eperp + math::cross(p_unit*
            static_cast<RealType>(light_speed), em_b));

        const auto coeff = static_cast<RealType>(1.0)/
            (static_cast<RealType>(electron_mass*light_speed*schwinger_field));

        return mod*norm_p*coeff;
    }

    //chi for photons
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_photon(
        const RealType px, const RealType py, const RealType pz,
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto p = math::vec3<RealType>{px, py, pz};
        const auto em_e = math::vec3<RealType>{ex, ey, ez};
        const auto em_b = math::vec3<RealType>{bx, by, bz};
        return chi_photon(p, em_e, em_b);
    }

    //chi for leptons
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_lepton(
        const math::vec3<RealType> t_p,
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType lambda = static_cast<RealType>(1.0))
    {

        const auto p = t_p * fact_momentum_to_SI_from<UnitSystem, RealType>();
        const auto em_e = t_em_e * fact_E_to_SI_from<UnitSystem, RealType>(lambda);
        const auto em_b = t_em_b * fact_B_to_SI_from<UnitSystem, RealType>(lambda);

        const auto norm_p = math::norm(p);
        if(norm_p == static_cast<RealType>(0.0))
            return static_cast<RealType>(0.0);

        const auto p_unit = p / norm_p;

        const auto one = static_cast<RealType>(1.0);
        const auto me_c = static_cast<RealType>(electron_mass*light_speed);

        //For gamma_2, writing the operations like this is better for single
        //precision (found with some tests).
        const auto norm_p_over_me_c = norm_p/me_c;
        const auto gamma_2 = one + (norm_p_over_me_c)*(norm_p_over_me_c);
        const auto gamma = sqrt(gamma_2);

        const auto beta = sqrt(one-one/gamma_2);
        const auto beta_vec = beta * p_unit;

        const auto beta_dot_e = math::dot(beta_vec, em_e);
        const auto beta_dot_e_2 = beta_dot_e*beta_dot_e;
        const auto e_plus_v_cross_b_2 = math::norm2(em_e + math::cross(
            beta_vec * static_cast<RealType>(light_speed), em_b));

        constexpr const auto coeff = static_cast<RealType>(1.0/schwinger_field);

        return gamma*sqrt(fabs(beta_dot_e_2-e_plus_v_cross_b_2))*coeff;
    }


    //chi for leptons
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_lepton(
        const RealType px, const RealType py, const RealType pz,
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto p = math::vec3<RealType>{px, py, pz};
        const auto em_e = math::vec3<RealType>{ex, ey, ez};
        const auto em_b = math::vec3<RealType>{bx, by, bz};
        return chi_photon(p, em_e, em_b);
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_CHI_FUNCTIONS
