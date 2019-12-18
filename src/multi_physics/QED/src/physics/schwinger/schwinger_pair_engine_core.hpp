#ifndef PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE

//This .hpp file contais the implementation of the
//Schwinger pair engine (as described in Gonoskov et al. PRE 92, 023305 2015
// and Banerjee et al. PRA 98, 032121 2018)

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

#include "schwinger_pair_engine_model.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{

    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int generate_schwinger_pairs_single(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        const RealType unf_zero_one_minus_epsi,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_schwinger_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto volume = dx*dy*dz;
        const auto probability = rate*volume*dt;

        return (unf_zero_one_minus_epsi < probability);
    }

    template<
        typename RealType,
        class RandomNumberGenerator,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int generate_schwinger_multiple_poisson(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        RandomNumberGenerator* rng,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_schwinger_pair_production_rate<RealType, UnitSystem>(
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
    int generate_schwinger_multiple_gaussian(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        RandomNumberGenerator* rng,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_schwinger_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto volume = dx*dy*dz;
        const auto probability = rate*volume*dt;

        const auto res = static_cast<RealType>(
            rng.gauss(probability, sqrt(probability));

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
    int generate_schwinger_multiple_choice(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType dx, const RealType dy, const RealType dz,
        const RealType dt,
        const RealType threshold,
        RandomNumberGenerator* rng,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto rate =
            compute_schwinger_pair_production_rate<RealType, UnitSystem>(
                ex, ey, ez, bx, by, bz, lambda);

        const auto volume = dx*dy*dz;
        const auto probability = rate*volume*dt;

        if(probability < static_cast<RealType>(threshold)){
            return rng->poisson(probability);
        }
        else{
            const auto res = static_cast<RealType>(
                rng.gauss(probability, sqrt(probability));

            if(res <= static_cast<RealType>(0.0))
                return 0;

            return static_cast<int>(res);
        }
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
