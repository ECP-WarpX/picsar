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
    * This function computes the number of expected Schwinger pairs
    * in a given volume for a given timestep using the Nikishov formula.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex the x component of the electric field
    * @param[in] ey the y component of the electric field
    * @param[in] ez the z component of the electric field
    * @param[in] bx the x component of the magnetic field
    * @param[in] by the y component of the magnetic field
    * @param[in] bz the z component of the magnetic field
    * @param[in] t_volume the volume of the physical region
    * @param[in] t_dt the temporal step
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the expected number of generated pairs
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_expected_pair_number(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType t_volume, const RealType t_dt,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        using namespace picsar::multi_physics::math;

        const auto ee = math::vec3<RealType>{ex,ey,ez}*
            fact_E_to_SI_from<UnitSystem, RealType>(ref_quantity);
        const auto bb = math::vec3<RealType>{bx,by,bz}*
            fact_B_to_SI_from<UnitSystem, RealType>(ref_quantity);

        constexpr const auto c_ff =
            static_cast<RealType>(1.0/(2.0*schwinger_field*schwinger_field));
        constexpr const auto c_gg =
            static_cast<RealType>(light_speed/schwinger_field/schwinger_field);

        const auto ff =
            (norm2(ee) - static_cast<RealType>(light_speed*light_speed)*norm2(bb))*c_ff;
        const auto gg = dot(ee, bb)*c_gg;

        const auto inner = sqrt(ff*ff+ gg*gg);

        const auto epsi = sqrt(inner + ff);
        const auto eta = sqrt(inner - ff);

        const auto zero = static_cast<RealType>(0.0);

        auto res = static_cast<RealType>(0.0);
        if(epsi != zero && eta != zero)
            res = epsi*eta*exp(-math::pi/epsi)/tanh(math::pi*eta/epsi);
        else if(epsi == zero)
            return zero;
        else
            res = epsi*epsi*exp(-math::pi/epsi)/math::pi;

        const auto dt = t_dt *
            fact_time_to_SI_from<UnitSystem>(ref_quantity);

        const auto volume = t_volume*
            fact_volume_to_SI_from<UnitSystem>(ref_quantity);

        constexpr const auto c1 = static_cast<RealType>(
            (elementary_charge*schwinger_field)*
            (elementary_charge*schwinger_field)/(4.0*math::pi*math::pi));

        const auto schwinger_pair_prod_coeff =
            static_cast<RealType>(c1*(volume/reduced_plank)*
            (dt/(reduced_plank*light_speed)));

        return schwinger_pair_prod_coeff*res;
}

    /**
    * This function computes the number of expected Schwinger pairs
    * and returns an integer number extracted from a Poisson distribution.
    * It requires a pointer to a RNG wrapper providing "poisson(lambda)" method
    * (see utils/rng_wrapper.hpp) for an example
    *
    * @tparam RealType the floating point type to be used
    * @tparam RandomNumberGenerator the RNG wrapper
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex the x component of the electric field
    * @param[in] ey the y component of the electric field
    * @param[in] ez the z component of the electric field
    * @param[in] bx the x component of the magnetic field
    * @param[in] by the y component of the magnetic field
    * @param[in] bz the z component of the magnetic field
    * @param[in] t_volume the volume of the physical region
    * @param[in] t_dt the temporal step
    * @param[in] unf_zero_one_minus_epsi a uniformly distributed random number
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the number of generated pairs
    */
    template<
        typename RealType,
        class RandomNumberGenerator,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int get_num_pairs_multiple_poisson(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType t_volume, const RealType t_dt,
        RandomNumberGenerator* rng,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        const auto expected_number =
            compute_expected_pair_number<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, t_volume, t_dt, ref_quantity);

        return rng->poisson(expected_number);
    }

    /**
    * This function computes the number of expected Schwinger pairs
    * and returns an integer number extracted from a Gaussian distribution.
    * It needs a normally distributed random number (mean=0, sigma=1)
    * Use this function only if the probability of generating
    * a pair is >> 1.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex the x component of the electric field
    * @param[in] ey the y component of the electric field
    * @param[in] ez the z component of the electric field
    * @param[in] bx the x component of the magnetic field
    * @param[in] by the y component of the magnetic field
    * @param[in] bz the z component of the magnetic field
    * @param[in] t_volume the volume of the physical region
    * @param[in] t_dt the temporal step
    * @param[in] gauss_mean_zero_sigma_one the normally distributed random number
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the number of generated pairs
    */
    template<
        typename RealType,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_num_pairs_multiple_gaussian(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType t_volume, const RealType t_dt,
        const RealType gauss_mean_zero_sigma_one,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        const auto expected_number =
            compute_expected_pair_number<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, t_volume, t_dt, ref_quantity);

        const auto num_pairs = expected_number +
            gauss_mean_zero_sigma_one*sqrt(expected_number);

        if(num_pairs <= static_cast<RealType>(0.0))
            return 0;

        return num_pairs;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
