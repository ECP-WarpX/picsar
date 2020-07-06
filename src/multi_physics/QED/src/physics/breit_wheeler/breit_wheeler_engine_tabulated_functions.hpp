#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS

//This .hpp file contais the implementation of the functions
//tabulated in the lookup tables for Breit-Wheeler pair production.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) A.I.Nikishov. & V.I. Ritus Sov. Phys. JETP 19, 2 (1964)
// 2) T.Erber Rev. Mod. Phys. 38, 626 (1966)
// 3) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 4) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses physical constants
#include "../phys_constants.h"
//Uses mathematical constants
#include "../../math/math_constants.h"
//Uses integration routines
#include "../../math/quadrature.hpp"
//Uses special functions
#include "../../math/special_functions.hpp"
//Uses sqrt and cbrt
#include "../../math/cmath_overloads.hpp"

#include <algorithm>
#include <cmath>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    /**
    * This function computes the X parameter (see validation script).
    * This function is not designed to be run on GPUs.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of the photon
    * @param[in] chi_ele the chi parameter of the electron
    *
    * @return the parameter X
    */
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType compute_x(
        const RealType chi_phot, const RealType chi_ele) noexcept
    {
        const auto temp = math::m_cbrt(chi_phot/(chi_ele*(chi_phot-chi_ele)));
        return temp*temp;
    }

    /**
    * This function computes the integrand for the
    * T function (see validation script). It is not usable
    * on GPUs.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of the photon
    * @param[in] chi_ele the chi parameter of the electron
    *
    * @return the value of the integrand for the T function
    */
    template<typename RealType>
    constexpr RealType compute_T_integrand(
        const RealType chi_phot, const RealType chi_ele)
    {
        using namespace math;
        if( chi_ele >= chi_phot )
            return zero<RealType>;

        if (chi_phot == zero<RealType> || chi_ele == zero<RealType>)
            return zero<RealType>;

        const auto xx = compute_x(chi_phot, chi_ele);
        const auto sqrt_xx = m_sqrt(xx);
        const auto xx_3_2 = sqrt_xx*sqrt_xx*sqrt_xx;

        const auto inner_integral = quad_a_inf<RealType>(
            [](RealType s){
                const auto sqrt_s = m_sqrt(s);
                const auto s_3_2 = sqrt_s*sqrt_s*sqrt_s;

                return sqrt_s*math::k_v(one_third<RealType>,
                    two_thirds<RealType>*s_3_2);
            }, xx);

        auto prod = (two<RealType>-chi_phot*xx_3_2)
            *k_v(two_thirds<RealType>,
                two_thirds<RealType>*xx_3_2);
        if(std::isnan(prod)) prod = zero<RealType>;

        return (inner_integral-prod);
    }

    /**
    * Computes the T function (see validation script).
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of the photon
    *
    * @return the value of the T function
    */
    template<typename RealType>
    RealType compute_T_function(const RealType chi_phot)
    {
        using namespace math;
        if(chi_phot <= math::zero<RealType>) return math::zero<RealType>;
        constexpr auto coeff = static_cast<RealType>(1./math::pi<>);
        return coeff*math::quad_a_b_s<RealType>(
            [=](RealType cc){
                return compute_T_integrand<RealType>(chi_phot, cc);},
            math::zero<RealType>, chi_phot)/(chi_phot*chi_phot*m_sqrt(3.0));
    }

    /**
    * Computes the numerator for the cumulative
    * probability (see validation script).
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of the photon
    * @param[in] chi_ele the chi parameter of the electron
    *
    * @return the value of the numerator of the cumulative probability distribution
    */
    template<typename RealType>
    RealType compute_cumulative_prob_numerator(
        const RealType chi_photon, RealType chi_ele)
    {
        using namespace math;
        constexpr auto coeff = static_cast<RealType>(1./pi<>);

        if(chi_photon <= zero<RealType>) return zero<RealType>;
        if(chi_ele  <= zero<RealType>) return zero<RealType>;

        if(chi_ele >= chi_photon) chi_ele = chi_photon;

        return coeff*quad_a_b_s<RealType>(
            [=](RealType cc){
                return compute_T_integrand<RealType>(chi_photon, cc);
            },zero<RealType>, chi_ele)/(chi_photon*chi_photon*m_sqrt(3.0));
    }

    /**
    * Computes the cumulative probability distribution
    * (see validation script) for a vector of chi paramters
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of the photon
    * @param[in] chis the chi parameters of the particles
    *
    * @return the cumulative probability distribution calculated for all the chi parameters
    */
    template<typename VectorType, typename RealType>
    VectorType compute_cumulative_prob(
        const RealType chi_photon, const VectorType& chis)
    {
        using namespace math;
        const auto den = compute_T_function(chi_photon);
        auto res = VectorType(chis.size());

        if(chi_photon <= zero<RealType>){
            for (auto& el: res ) el = zero<RealType>;
            return res;
        }

        //This may happen for chi_photon very small and single precision
        if (den <= zero<RealType>){
            std::transform(chis.begin(), chis.end(),
            res.begin(), [=](auto chi_part){
                if(chi_part < half<RealType>*chi_photon - std::numeric_limits<RealType>::epsilon())
                    return zero<RealType>;
                if(chi_part > half<RealType>*chi_photon + std::numeric_limits<RealType>::epsilon())
                    return one<RealType>;
                return half<RealType>;
            });
            return res;
        }

        std::transform(chis.begin(), chis.end(),
            res.begin(), [=](auto chi_part){
                const auto val =
                    compute_cumulative_prob_numerator(chi_photon, chi_part)/den;
                if(val <= one<RealType>) return val;
                return one<RealType>;});
        return res;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABULATED_FUNCTIONS
