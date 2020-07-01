#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS

#include <cmath>
#include <algorithm>

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
        if(chi_phot <= math::zero<RealType>) return math::zero<RealType>;
        constexpr auto coeff = static_cast<RealType>(1./math::pi<>);
        return coeff*math::quad_a_b<RealType>(
            [=](RealType cc){
                return compute_T_integrand<RealType>(chi_phot, cc);},
            math::zero<RealType>, chi_phot)/(chi_phot*chi_phot*sqrt(3.0));
    }

    template<typename RealType>
    RealType compute_cumulative_prob_numerator(
        const RealType chi_photon, RealType chi_part)
    {
        constexpr auto coeff = static_cast<RealType>(1./math::pi<>);

        if(chi_photon <= math::zero<RealType>) return math::zero<RealType>;
        if(chi_part  <= math::zero<RealType>) return math::zero<RealType>;

        if(chi_part >= chi_photon) chi_part = chi_photon;

        return coeff*math::quad_a_b<RealType>(
            [=](RealType cc){
                return compute_T_integrand<RealType>(chi_photon, cc);
            },math::zero<RealType>, chi_part)/(chi_photon*chi_photon*sqrt(3.0));
    }

    template<typename VectorType, typename RealType>
    VectorType compute_cumulative_prob(
        const RealType chi_photon, const VectorType& chi_particle)
    {
        const auto den = compute_T_function(chi_photon);
        auto res = VectorType(chi_particle.size());

        if(chi_photon <= math::zero<RealType>){
            for (auto& el: res ) el = math::zero<RealType>;
            return res;
        }

        //this may happen for chi_photon very small and single precision
        if (den <= math::zero<RealType>){
            std::transform(chi_particle.begin(), chi_particle.end(),
            res.begin(), [=](auto chi_part){
                if(chi_part < math::half<RealType>*chi_photon - std::numeric_limits<RealType>::epsilon())
                    return math::zero<RealType>;
                if(chi_part > math::half<RealType>*chi_photon + std::numeric_limits<RealType>::epsilon())
                    return math::one<RealType>;
                return math::half<RealType>;
            });
            return res;
        }

        std::transform(chi_particle.begin(), chi_particle.end(),
            res.begin(), [=](auto chi_part){
                const auto val =
                    compute_cumulative_prob_numerator(chi_photon, chi_part)/den;
                if(val <= math::one<RealType>) return val;
                return math::one<RealType>;});
        return res;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
