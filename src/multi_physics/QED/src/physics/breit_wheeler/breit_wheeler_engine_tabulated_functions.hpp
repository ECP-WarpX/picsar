#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS

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
namespace breit_wheeler{

    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType compute_x(
        const RealType chi_phot, const RealType chi_ele) noexcept
    {
        const auto temp = cbrt(chi_phot/(chi_ele*(chi_phot-chi_ele)));
        return temp*temp;
    }

    template<typename RealType>
    constexpr RealType compute_T_integrand(
        const RealType chi_phot, const RealType chi_ele) noexcept
    {
        using namespace math;
        if( chi_ele >= chi_phot )
            return zero<RealType>;

        if (chi_phot == zero<RealType> || chi_ele == zero<RealType>)
            return zero<RealType>;

        const auto xx = compute_x(chi_phot, chi_ele);
        const auto sqrt_xx = sqrt(xx);
        const auto xx_3_2 = sqrt_xx*sqrt_xx*sqrt_xx;

        const auto inner_integral = math::quad_a_inf<RealType>(
            [](RealType s){
                const auto sqrt_s = sqrt(s);
                const auto s_3_2 = sqrt_s*sqrt_s*sqrt_s;

                return sqrt_s*math::k_v(math::one_third<RealType>,
                    math::two_thirds<RealType>*s_3_2);
            }, xx);

        return (inner_integral-(math::two<RealType>-chi_phot*xx_3_2)
            *k_v(math::two_thirds<RealType>,
                math::two_thirds<RealType>*xx_3_2));
    }

    template<typename RealType>
    RealType compute_T_function(const RealType chi_phot)
    {
        constexpr auto coeff = static_cast<RealType>(1./math::pi<>);
        return coeff*math::quad_a_b<RealType>(
            [=](RealType cc){
                return compute_T_integrand<RealType>(chi_phot, cc);},
            math::zero<RealType>, chi_phot)/(chi_phot*chi_phot*sqrt(3.0));
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
