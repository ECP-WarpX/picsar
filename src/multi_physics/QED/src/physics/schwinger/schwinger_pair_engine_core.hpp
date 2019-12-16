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

//############################################### Declaration

namespace picsar{
namespace multi_physics{
namespace phys{

    //Computes the pair production rate per unit time per unit volume
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_schwinger_pair_production_rate(
        const RealType ex, const RealType ey, const RealType ez,
        const RealType bx, const RealType by, const RealType bz,
        const RealType lambda = static_cast<RealType>(1.0))
    {
        const auto ee = math::vec3<RealType>{ex,ey,ez}*
            fact_E_to_SI_from<UnitSystem, RealType>(lambda);
        const auto bb = math::vec3<RealType>{bx,by,bz}*
            fact_B_to_SI_from<UnitSystem, RealType>(lambda);

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
            fact_length_from_SI_to<UnitSystem, RealType>(lambda),3));

        constexpr const auto coeff = schwinger_pair_prod_coeff*
            fact_rate_from_SI_to<UnitSystem, RealType>(lambda)/unit_vol;

        return coeff*res;
}

    //This function determines if a pair has been generated in a given
    //cell and its statistical weight.
    //It returns a bool (true if a pair has been generated)
    //and a _REAL (the statistical weight of the new pair)
    //It requires to provide the fields, the cell size and dt
    //Use this function if the probability to generate a pair is << 1
    template<typename RealType, unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    int generate_pairs_single(
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

}
}
}



//______________________GPU
//Same as above, but with GPU use directly this one!
//Returns false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
void
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
internal_generate_pairs_multiple(
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dx, _REAL dy, _REAL dz,
_REAL dt,
size_t* how_many, _REAL* weight,
_REAL _lambda, _REAL unf_zero_one_minus_epsi
)
{
#ifdef PXRMP_WITH_SI_UNITS
    _lambda = static_cast<_REAL>(1.0);
#endif

    _REAL rate = internal_compute_schwinger_pair_production_rate
        (ex, ey, ez, bx, by, bz, _lambda);

    _REAL volume = dx*dy*dz;
    _REAL probability = rate*volume*dt;

    *weight = one / volume;

    *how_many = poisson_distrib(probability, unf_zero_one_minus_epsi);
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
