#ifndef PICSAR_MULTIPHYSICS_CHI_FUNCTIONS
#define PICSAR_MULTIPHYSICS_CHI_FUNCTIONS

//Should be included by all the src files of the library
#include "../qed_commons.h"

//Uses GPU-friendly array
#include "../containers/picsar_array.hpp"
//Uses operations on 3-vectors
#include "../math/vec_functions.hpp"
//Uses unit conversion
#include "unit_conversion.hpp"
//Uses sqrt
#include "../math/cmath_overloads.hpp"

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
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = t_p *conv<
            quantity::momentum, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact();
        const auto em_e = t_em_e * conv<
            quantity::E, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(reference_quantity);
        const auto em_b = t_em_b * conv<
            quantity::B, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(reference_quantity);

        const auto norm_p = norm(p);
        if(norm_p == zero<RealType>)
            return zero<RealType>;

        const auto p_unit = p / norm_p;
        const auto em_eperp = em_e - dot(p_unit,em_e)*p_unit;
        const auto field = norm(em_eperp + cross(p_unit, em_b));

        const auto gamma_phot = norm_p/
            heaviside_lorentz_electron_rest_energy<RealType>;

        constexpr auto one_over_schwinger = static_cast<RealType>(1.0/
            heaviside_lorentz_schwinger_field<double>);

        return gamma_phot*field*one_over_schwinger;
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
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = vec3<RealType>{px, py, pz};
        const auto em_e = vec3<RealType>{ex, ey, ez};
        const auto em_b = vec3<RealType>{bx, by, bz};
        return chi_photon<RealType, UnitSystem>(
            p, em_e, em_b, reference_quantity);
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
    RealType chi_ele_pos(
        const math::vec3<RealType> t_p,
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = t_p * conv<
            quantity::momentum, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact();
        const auto em_e = t_em_e * conv<
            quantity::E, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(reference_quantity);
        const auto em_b = t_em_b * conv<
            quantity::B, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(reference_quantity);

        const auto norm_p = norm(p);
        if(norm_p == zero<RealType>)
            return zero<RealType>;

        const auto p_unit = p / norm_p;

        constexpr auto one_over_m2 = static_cast<RealType>(1.0/
            heaviside_lorentz_electron_rest_energy<double>/
            heaviside_lorentz_electron_rest_energy<double>);
        const auto gamma_2 = one<RealType> + norm_p*norm_p*one_over_m2;
        const auto gamma = m_sqrt(gamma_2);

        const auto beta =m_sqrt((gamma_2-one<RealType>)/gamma_2);
        const auto beta_vec = beta * p_unit;

        const auto beta_dot_e = dot(beta_vec, em_e);
        const auto beta_dot_e_2 = beta_dot_e*beta_dot_e;
        const auto e_plus_beta_cross_b_2 =
            norm2(em_e + cross(beta_vec, em_b));
        const auto field = m_sqrt(fabs(beta_dot_e_2-e_plus_beta_cross_b_2));

        constexpr auto one_over_schwinger = static_cast<RealType>(1.0/
            heaviside_lorentz_schwinger_field<double>);

        return gamma*field*one_over_schwinger;
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
    RealType chi_ele_pos(
        const RealType px, const RealType py, const RealType pz,
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType reference_quantity = math::one<RealType>)
    {
        const auto p = math::vec3<RealType>{px, py, pz};
        const auto em_e = math::vec3<RealType>{ex, ey, ez};
        const auto em_b = math::vec3<RealType>{bx, by, bz};
        return chi_ele_pos<RealType, UnitSystem>(
            p, em_e, em_b, reference_quantity);
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_CHI_FUNCTIONS
