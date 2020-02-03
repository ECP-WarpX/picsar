#ifndef PICSAR_MULTIPHYSICS_CHI_FUNCTIONS
#define PICSAR_MULTIPHYSICS_CHI_FUNCTIONS

#include <cmath>

//Should be included by all the src files of the library
#include "../qed_commons.h"

//Uses picsar_arrays
#include "../containers/picsar_array.hpp"

#include "../math/vec_functions.hpp"

#include "unit_conversion.hpp"

// This file contains auxiliary functions to calculate the chi parameter
// for photons and leptons (electrons and positrons)

namespace picsar{
namespace multi_physics{
namespace phys{

    /**
    * This function returns the chi parameter for a photon
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input&output [default is SI]
    * @param[in] t_p a vec3<RealType> containing the momentum of the photon
    * @param[in] t_em_e a vec3<RealType> containing the electric field
    * @param[in] t_em_b a vec3<RealType> containing the magnetic field
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the chi parameter
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_photon(
        const math::vec3<RealType> t_p,
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {

        const auto p = t_p * fact_momentum_to_SI_from<UnitSystem, RealType>();
        const auto em_e = t_em_e *
            fact_E_to_SI_from<UnitSystem, RealType>(reference_quantity);
        const auto em_b = t_em_b *
            fact_B_to_SI_from<UnitSystem, RealType>(reference_quantity);

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

    /**
    * This function returns the chi parameter for a photon
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input&output [default is SI]
    * @param[in] px a RealType containing the x component of the momentum of the photon
    * @param[in] py a RealType containing the y component of the momentum of the photon
    * @param[in] pz a RealType containing the z component of the momentum of the photon
    * @param[in] ex a RealType containing the x component of the electric field
    * @param[in] ey a RealType containing the y component of the electric field
    * @param[in] ez a RealType containing the z component of the electric field
    * @param[in] bx a RealType containing the x component of the magnetic field
    * @param[in] by a RealType containing the y component of the magnetic field
    * @param[in] bz a RealType containing the z component of the magnetic field
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the chi parameter
    */
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_photon(
        const RealType px, const RealType py, const RealType pz,
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        const auto p = math::vec3<RealType>{px, py, pz};
        const auto em_e = math::vec3<RealType>{ex, ey, ez};
        const auto em_b = math::vec3<RealType>{bx, by, bz};
        return chi_photon(p, em_e, em_b, reference_quantity);
    }

    /**
    * This function returns the chi parameter for either an electron or a positron
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input&output [default is SI]
    * @param[in] t_p a vec3<RealType> containing the momentum of the electron or the positrion
    * @param[in] t_em_e a vec3<RealType> containing the electric field
    * @param[in] t_em_b a vec3<RealType> containing the magnetic field
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the chi parameter
    */
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_lepton(
        const math::vec3<RealType> t_p,
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {

        const auto p = t_p * fact_momentum_to_SI_from<UnitSystem, RealType>();
        const auto em_e = t_em_e * fact_E_to_SI_from<UnitSystem, RealType>(reference_quantity);
        const auto em_b = t_em_b * fact_B_to_SI_from<UnitSystem, RealType>(reference_quantity);

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

        const auto beta = static_cast<RealType>(sqrt((gamma_2-one)/gamma_2));
        const auto beta_vec = beta * p_unit;

        const auto beta_dot_e = math::dot(beta_vec, em_e);
        const auto beta_dot_e_2 = beta_dot_e*beta_dot_e;
        const auto e_plus_v_cross_b_2 = math::norm2(em_e + math::cross(
            beta_vec * static_cast<RealType>(light_speed), em_b));

        constexpr const auto coeff = static_cast<RealType>(1.0/schwinger_field);

        return gamma*sqrt(fabs(beta_dot_e_2-e_plus_v_cross_b_2))*coeff;
    }


    /**
    * This function returns the chi parameter for either an electron or a positron
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input&output [default is SI]
    * @param[in] px a RealType containing the x component of the momentum of the electron (or the positron)
    * @param[in] py a RealType containing the y component of the momentum of the electron (or the positron)
    * @param[in] pz a RealType containing the z component of the momentum of the electron (or the positron)
    * @param[in] ex a RealType containing the x component of the electric field
    * @param[in] ey a RealType containing the y component of the electric field
    * @param[in] ez a RealType containing the z component of the electric field
    * @param[in] bx a RealType containing the x component of the magnetic field
    * @param[in] by a RealType containing the y component of the magnetic field
    * @param[in] bz a RealType containing the z component of the magnetic field
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the chi parameter
    */
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType chi_lepton(
        const RealType px, const RealType py, const RealType pz,
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        const auto p = math::vec3<RealType>{px, py, pz};
        const auto em_e = math::vec3<RealType>{ex, ey, ez};
        const auto em_b = math::vec3<RealType>{bx, by, bz};
        return chi_photon(p, em_e, em_b, reference_quantity);
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_CHI_FUNCTIONS
