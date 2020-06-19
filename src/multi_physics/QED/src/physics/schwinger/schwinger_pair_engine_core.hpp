#ifndef PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE

//This .hpp file contais the implementation of the
//core functions of the Schwinger pair engine.
//If desired, these functions can be used directly
//without the higher-level interface.
//
// References:
// 1) Schwinger. Phys. Rev. 82, 5 (1951)
// 2) Nikishov. Sov. Phys. JETP 30, 660 (1970)
// 3) Narozhny et al. Phys. LEtt. A, 330, 1-2 (2004)
// 4) Bulanov et al. Phys. Rev. Lett. 104, 220404 (2010)

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar arrays
#include "../../math/vec_functions.hpp"

//Uses vector functions
#include "../../math/vec_functions.hpp"

//Uses physical constants
#include "../phys_constants.h"

//Uses unit conversion"
#include "../unit_conversion.hpp"

//Uses Poisson's distribution
#include "../../utils/rng_wrapper.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace schwinger{

    /**
    * This function computes the Schwinger pair production rate
    * using the Nikishov formula.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] t_em_e a vec3<RealType> containing the electric field
    * @param[in] t_em_b a vec3<RealType> containing the magnetic field
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the pair production rate
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType pair_production_rate(
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        using namespace picsar::multi_physics::math;

        const auto em_e = t_em_e * conv<
            quantity::E, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);
        const auto em_b = t_em_b * conv<
            quantity::B, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        constexpr auto one_half = static_cast<RealType>(0.5);
        constexpr auto one_over_schwinger = static_cast<RealType>(
            1.0/heaviside_lorentz_schwinger_field<double>);

        const auto ff =(norm2(em_e) - norm2(em_b))*one_half;
        const auto gg = dot(em_e, em_b);

        const auto inner = sqrt(ff*ff+ gg*gg);

        const auto epsi = sqrt(fabs(inner + ff))*one_over_schwinger;
        const auto eta = sqrt(fabs(inner - ff))*one_over_schwinger;

        constexpr const auto coeff = static_cast<RealType>(
            heaviside_lorentz_elementary_charge<double>*
            heaviside_lorentz_elementary_charge<double>*
            heaviside_lorentz_schwinger_field<double>*
            heaviside_lorentz_schwinger_field<double>/
            (4.0*pi<double>*pi<double>));
        const auto rate_conv = conv<quantity::rate,
            unit_system::heaviside_lorentz, UnitSystem, RealType>::fact(1.0, ref_quantity);

        if(epsi != zero<RealType> && eta != zero<RealType>)
            return coeff*rate_conv*epsi*eta*exp(-pi<RealType>/epsi)/tanh(pi<RealType>*eta/epsi);
        else if(epsi == zero<RealType>)
            return zero<RealType>;
        else
            return coeff*rate_conv*epsi*epsi*exp(-pi<RealType>/epsi)/pi<RealType>;
    }

    /**
    * This function computes the Schwinger pair production rate
    * using the Nikishov formula.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex a RealType containing the x component of the electric field
    * @param[in] ey a RealType containing the y component of the electric field
    * @param[in] ez a RealType containing the z component of the electric field
    * @param[in] bx a RealType containing the x component of the magnetic field
    * @param[in] by a RealType containing the y component of the magnetic field
    * @param[in] bz a RealType containing the z component of the magnetic field
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the pair production rate
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType pair_production_rate(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType ref_quantity = math::one<RealType>)
    {
        using namespace math;
        const auto em_e = vec3<RealType>{ex, ey, ez};
        const auto em_b = vec3<RealType>{bx, by, bz};
        return pair_production_rate<RealType, UnitSystem>(em_e, em_b,
            ref_quantity);
    }

    /**
    * This function computes the number of expected Schwinger pairs
    * in a given volume for a given timestep using the Nikishov formula.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] t_em_e a vec3<RealType> containing the electric field
    * @param[in] t_em_b a vec3<RealType> containing the magnetic field
    * @param[in] t_volume the volume of the physical region
    * @param[in] t_dt the temporal step
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the expected number of generated pairs
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType expected_pair_number(
        const math::vec3<RealType> t_em_e,
        const math::vec3<RealType> t_em_b,
        const RealType t_volume, const RealType t_dt,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        using namespace picsar::multi_physics::math;

        const auto em_e = t_em_e * conv<
            quantity::E, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);
        const auto em_b = t_em_b * conv<
            quantity::B, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        const auto schwinger_pair_production_rate = pair_production_rate<
            RealType, unit_system::heaviside_lorentz>(em_e, em_b);

        const auto dt = t_dt * conv<
            quantity::time, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);
        const auto volume = t_volume * conv<
            quantity::volume, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        return schwinger_pair_production_rate*dt*volume;
    }

    /**
    * This function computes the number of expected Schwinger pairs
    * in a given volume for a given timestep using the Nikishov formula.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex a RealType containing the x component of the electric field
    * @param[in] ey a RealType containing the y component of the electric field
    * @param[in] ez a RealType containing the z component of the electric field
    * @param[in] bx a RealType containing the x component of the magnetic field
    * @param[in] by a RealType containing the y component of the magnetic field
    * @param[in] bz a RealType containing the z component of the magnetic field
    * @param[in] t_volume the volume of the physical region
    * @param[in] t_dt the temporal step
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the expected number of generated pairs
    */
    template<typename RealType, unit_system UnitSystem>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType expected_pair_number(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType t_volume, const RealType t_dt,
        const RealType reference_quantity = math::one<RealType>)
    {
        const auto em_e = math::vec3<RealType>{ex, ey, ez};
        const auto em_b = math::vec3<RealType>{bx, by, bz};
        return expected_pair_number<RealType, UnitSystem>(
            em_e, em_b, t_volume, t_dt, reference_quantity);
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
