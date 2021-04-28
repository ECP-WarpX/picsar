#ifndef PICSAR_MULTIPHYSICS_GAMMA_FUNCTIONS
#define PICSAR_MULTIPHYSICS_GAMMA_FUNCTIONS

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Uses GPU-friendly array
#include "picsar_qed/containers/picsar_array.hpp"
//Uses operations on 3-vectors
#include "picsar_qed/math/vec_functions.hpp"
//Uses unit conversion
#include "picsar_qed/physics/unit_conversion.hpp"
//Uses sqrt
#include "picsar_qed/math/cmath_overloads.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{

    /**
    * This function returns the energy of a photon normalized with respect to the rest mass of an electron
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input [default is SI]
    * @param[in] t_p a vec3<RealType> containing the momentum of the photon
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the normalized energy of the photon
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    RealType compute_gamma_photon(
        const math::vec3<RealType> t_p,
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = t_p *conv<
            quantity::momentum, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(reference_quantity);

        const auto norm_p = norm(p);
        if(norm_p == zero<RealType>)
            return zero<RealType>;

        const auto gamma_phot = norm_p/
            heaviside_lorentz_electron_rest_energy<RealType>;

        return gamma_phot;
    }

    /**
    * This function returns the energy of a photon normalized with respect to the rest mass of an electron
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input [default is SI]
    * @param[in] px a RealType containing the x component of the momentum of the photon
    * @param[in] py a RealType containing the y component of the momentum of the photon
    * @param[in] pz a RealType containing the z component of the momentum of the photon
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the normalized energy of the photon
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    RealType compute_gamma_photon(
        const RealType px, const RealType py, const RealType pz,
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = vec3<RealType>{px, py, pz};
        return compute_gamma_photon<RealType, UnitSystem>(p, reference_quantity);
    }

    /**
    * This function returns the Lorentz factor for an electron or a positron
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input [default is SI]
    * @param[in] t_p a vec3<RealType> containing the momentum of the electron (or positron)
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the Lorentz factor
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    RealType compute_gamma_ele_pos(
        const math::vec3<RealType> t_p,
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = t_p *conv<
            quantity::momentum, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(reference_quantity);

        const auto norm_p = norm(p);

        constexpr auto one_over_m2 = static_cast<RealType>(1.0/
            heaviside_lorentz_electron_rest_energy<double>/
            heaviside_lorentz_electron_rest_energy<double>);
        const auto gamma_2 = one<RealType> + norm_p*norm_p*one_over_m2;
        const auto gamma = m_sqrt(gamma_2);

        return gamma;
    }

    /**
    * This function returns the Lorentz factor for an electron or a positron
    *
    * @tparam RealType the floating point type of the result
    * @tparam unit_system unit system to be used for input [default is SI]
    * @param[in] px a RealType containing the x component of the momentum of the electron (or the positron)
    * @param[in] py a RealType containing the y component of the momentum of the electron (or the positron)
    * @param[in] pz a RealType containing the z component of the momentum of the electron (or the positron)
    * @param[in] reference_quantity reference quantity for unit conversion (e.g. lambda or omega_r), if needed
    * @return the Lorentz factor
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    RealType compute_gamma_ele_pos(
        const RealType px, const RealType py, const RealType pz,
        const RealType reference_quantity = math::one<RealType>)
    {
        using namespace math;

        const auto p = math::vec3<RealType>{px, py, pz};
        return compute_gamma_ele_pos<RealType, UnitSystem>(p, reference_quantity);
    }

}
}
}

#endif //PICSAR_MULTIPHYSICS_GAMMA_FUNCTIONS
