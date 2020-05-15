#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS

#include <cmath>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses physical constants
#include "../phys_constants.h"

//Uses physical constants
#include "../../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType compute_x(
        const RealType chi_phot, const RealType chi_ele) noexcept
    {
        const auto temp = cbrt(chi_phot/(chi_ele*(chi_phot-chi_ele));
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

        constexpr auto div = one<RealType>/(pi<RealType>*three<RealType>);


        constexpr auto one_third = one<RealType>/three<RealType>;
        constexpr auto two_thirds = one<RealType>/three<RealType>;
        const auto inner_integral = math::quad_a_inf<RealType>(
            [](RealType s){
                const auto sqrt_s = sqrt(s);
                const auto s_3_2 = sqrt_s*sqrt_s*sqrt_s;

                return sqrt_s*k_v(one_third<RealType>, two_thirds<RealType>*s_3_2);
            }, xx);

        return (inner_integral-(two<RealType>-chi_phot*xx_3_2)
            *k_v(two_thirds<RealType>, two_thirds<RealType>*xx_3_2))*div;
    }

    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType compute_T_function(const RealType chi_phot)
    {
        using namespace math;
        return quad_a_b<RealType>(compute_T_integrand<RealType>, zero<RealType>, chi_phot);
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
