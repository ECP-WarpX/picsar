#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABULATED_FUNCTIONS
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABULATED_FUNCTIONS

#include <cmath>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses physical constants
#include "../phys_constants.h"

//Uses physical constants
#include "../../math/math_constants.h"

#include "../../math/quadrature.hpp"

#include "../../math/special_functions.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType compute_y(
        const RealType chi_part, const RealType csi) noexcept
    {

        return math::two_thirds<RealType>*
            csi/(chi_part*(math::one<RealType> - csi));
    }

    template<typename RealType>
    constexpr RealType compute_G_integrand(
        const RealType chi_part, const RealType chi_phot) noexcept
    {
        using namespace math;
        if( chi_phot >= chi_part )
            return zero<RealType>;

        if (chi_phot == zero<RealType> || chi_part == zero<RealType>)
            return zero<RealType>;

        const auto csi = chi_phot/chi_part;

        const auto yy = compute_y(chi_part, csi);

        constexpr RealType coeff = sqrt(three<>)/(two<>*pi<>);

        const auto inner_integral = quad_a_inf<RealType>(
            [](RealType s){return k_v(five_thirds<RealType>,s);}, yy);

        const auto second_part = (csi*csi/(one<RealType>-csi))*
            k_v(two_thirds<RealType>,yy);

        return coeff*(inner_integral + second_part)/chi_part;
    }

    template<typename RealType>
    RealType compute_G_function(const RealType chi_part,
        const RealType chi_phot_min = math::zero<RealType>)
    {
        using namespace math;
        return quad_a_b<RealType>(
            [=](RealType chi_phot){
                return compute_G_integrand<RealType>(chi_part, chi_phot);},
                chi_phot_min, chi_part);
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABULATED_FUNCTIONS
