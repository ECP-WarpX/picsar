#ifndef PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE

//This .hpp file contais the implementation of the
//core functions of the Schwinger pair engine.
//If desired, these functions can be used directly
//without the higher-level interface.
//
// References:
// 1) Gonoskov et al. PRE 92, 023305 2015
// 2) Banerjee et al. PRA 98, 032121 2018

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

namespace picsar{
namespace multi_physics{
namespace phys{
namespace schwinger{

    /**
    * This function computes the Schwinger pair production
    * rate per unit time per unit volume.
    *
    * @tparam RealType the floating point type to be used
    * @tparam UnitSystem unit system to be used for inputs & outputs
    * @param[in] ex the x component of the electric field
    * @param[in] ey the y component of the electric field
    * @param[in] ez the z component of the electric field
    * @param[in] bx the x component of the magnetic field
    * @param[in] by the y component of the magnetic field
    * @param[in] bz the z component of the magnetic field
    * @param[in] ref_quantity reference quantityt for unit conversion (lambda or omega)

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

        constexpr const auto unit_vol = static_cast<RealType>(pow(
            fact_length_from_SI_to<UnitSystem, RealType>(ref_quantity),3));

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
    * This function computes the Schwinger pair production
    * and returns 1 if a pair is generated, 0 otherwise.
    * It needs a uniformly distributed random number.
    * Use this function only if the probability of generating
    * a pair is << 1.
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
    * @param[in] dt the temporal step
    * @param[in] unf_zero_one_minus_epsi a uniformly distributed random number
    * @param[in] ref_quantity reference quantityt for unit conversion (lambda or omega)

    * @return 1 if a particle is generated, 0 otherwise.
    */
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int get_num_pairs_single(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType t_volume, const RealType dt,
        const RealType unf_zero_one_minus_epsi,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, ref_quantity);

        const auto volume = t_volume;
        const auto probability = rate*volume*dt;

        return (unf_zero_one_minus_epsi < probability);
    }

    template<
        typename RealType,
        class RandomNumberGenerator,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int get_num_pairs_multiple_poisson(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        RandomNumberGenerator* rng,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto volume = dx*dy*dz;
        const auto probability = rate*volume*dt;

        return rng->poisson(probability);
    }

    template<
        typename RealType,
        class RandomNumberGenerator,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int get_num_pairs_multiple_gaussian(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        RandomNumberGenerator* rng,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto volume = dx*dy*dz;
        const auto probability = rate*volume*dt;

        const auto res = static_cast<RealType>(
            rng->gaussian(probability, sqrt(probability)));

        if(res <= static_cast<RealType>(0.0))
            return 0;

        return static_cast<int>(res);
    }

    //This function determines how many pairs have been generated in a given
    //cell
    template<
        typename RealType,
        class RandomNumberGenerator,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int get_num_pairs_multiple_choice(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        const RealType threshold,
        RandomNumberGenerator* rng,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto volume = dx*dy*dz;
        const auto probability = rate*volume*dt;

        if(probability < static_cast<RealType>(threshold)){
            return rng->poisson(probability);
        }
        else{
            const auto res = static_cast<RealType>(
                rng->gaussian(probability, sqrt(probability)));

            if(res <= static_cast<RealType>(0.0))
                return 0;

            return static_cast<int>(res);
        }
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
