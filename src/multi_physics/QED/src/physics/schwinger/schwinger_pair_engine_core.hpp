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
// 5) Gonoskov et al. Phys. Rev. E 92, 023305 (2015)


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
    * This function computes the Schwinger pair production
    * rate per unit time per unit volume using the Nikishov formula.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex the x component of the electric field
    * @param[in] ey the y component of the electric field
    * @param[in] ez the z component of the electric field
    * @param[in] bx the x component of the magnetic field
    * @param[in] by the y component of the magnetic field
    * @param[in] bz the z component of the magnetic field
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the conversion factor
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_pair_production_rate(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
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

        constexpr const auto unit_vol =
            fact_volume_from_SI_to<UnitSystem, RealType>(ref_quantity);

        constexpr const auto schwinger_pair_prod_coeff =
            static_cast<RealType>{
                elementary_charge*elementary_charge*
                schwinger_field*schwinger_field/
                (4.0*math::pi*math::pi*reduced_plank*reduced_plank*light_speed)};

        constexpr const auto coeff = schwinger_pair_prod_coeff*
            fact_rate_from_SI_to<UnitSystem, RealType>(ref_quantity)/unit_vol;

        return coeff*res;
}

    /**
    * This function computes the Schwinger pair production rate
    * and returns an integer number extracted from a Poisson's distribution
    * according to the calculated rate.
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
        const auto dt = t_dt *
            fact_time_to_SI_from<UnitSystem>(ref_quantity);

        const auto volume = t_volume*
            fact_volume_to_SI_from<UnitSystem>(ref_quantity);

        const auto rate =
            compute_pair_production_rate<RealType, unit_system::SI>(
                ex, ey, ez, bx, by, bz);

        const auto probability = rate*volume*dt;

        return rng->poisson(probability);
    }

    /**
    * This function computes the Schwinger pair production rate
    * and returns a real number extracted from a Gaussian distribution
    * according to the calculated rate.
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
        const RealType t_volume, const RealType dt,
        const RealType gauss_mean_zero_sigma_one,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        const auto dt = t_dt *
            fact_time_to_SI_from<UnitSystem>(ref_quantity);

        const auto volume = t_volume*
            fact_volume_to_SI_from<UnitSystem>(ref_quantity);

        const auto rate =
            compute_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto probability = rate*volume*dt;

        const auto num_pairs = probability +
            gauss_mean_zero_sigma_one*sqrt(probability);

        if(num_pairs <= static_cast<RealType>(0.0))
            return 0;

        return num_pairs;
    }

    /**
    * This function computes the Schwinger pair production rate
    * and returns a real number extracted from a Gaussian distribution
    * according to the calculated rate. Indeed, if the pair prodduction rate
    * is >> 1, a Gaussian approximation can be used. A threshold
    * to choose between the two distributions can be specified by the user
    * (the default value is 1000). This function needs a RNG wrapper able
    * to provide the functions  "poisson(lambda)" and "gaussian(mean, deviation)".
    * See utils/rng_wrapper.hpp for an example.
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
    * @param[in] rng a pointer to the RNG wrapper
    * @param[in] ref_quantity reference quantity for unit conversion (lambda or omega)

    * @return the number of generated pairs
    */
    template<
        typename RealType,
        class RandomNumberGenerator,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_num_pairs_multiple_choice(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType t_volume, const RealType dt,
        RandomNumberGenerator* rng,
        const RealType threshold = static_cast<RealType>(1000.),
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        const auto dt = t_dt *
            fact_time_to_SI_from<UnitSystem>(ref_quantity);

        const auto volume = t_volume*
            fact_volume_to_SI_from<UnitSystem>(ref_quantity);

        const auto rate =
            compute_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto probability = rate*volume*dt;

        if(probability < static_cast<RealType>(threshold)){
            return static_cast<RealType>(rng->poisson(probability);
        }
        else{
            const auto res = static_cast<RealType>(
                rng->gaussian(probability, sqrt(probability)));

            if(res <= static_cast<RealType>(0.0))
                return 0;

            return res;
        }
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
