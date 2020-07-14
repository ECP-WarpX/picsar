#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABULATED_FUNCTIONS
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABULATED_FUNCTIONS

//This .hpp file contais the implementation of the functions
//tabulated in the lookup tables for the Quantum Synchrotron photon emission engine.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

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

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    /**
    * It computes the y parameter (see validation script).
    * This function is not designed to be run on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_phot the chi parameter of the photon
    * @param[in] chi_ele the chi parameter of the electron
    *
    * @return the parameter X
    */
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType compute_y(
        const RealType chi_part, const RealType csi) noexcept
    {
        return math::two_thirds<RealType>*
            csi/(chi_part*(math::one<RealType> - csi));
    }

    /**
    * This function replaces the difficult integration of
    * kv(5/3, x) from y to infinity with a rapidly convergening expression
    * (see https://en.wikipedia.org/wiki/Synchrotron_radiation
    * and M.Kh.Khokonov. Journal of Experimental and Theoretical
    * Physics volume 99, 690–707(2004))
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] y the left extremum of the integration
    *
    * @return the integral of kv(5/3, x) from y to infinity
    */
    template<typename RealType>
    RealType inner_integral(const RealType y)
    {
        using namespace math;
        return quad_a_inf<RealType>(
            [=](RealType s){
                using namespace math;
                if( y ==  zero<RealType>)
                    return std::numeric_limits<RealType>::infinity();;
                const auto s2 = s*s;
                const auto s4 = s2*s2;
                const auto cc = (one<RealType> +
                    four<RealType>*one_third<RealType>*s2)*
                    m_sqrt(one<RealType> + one_third<RealType>*s2);
                const auto f1 = static_cast<RealType>(9.0) +
                    static_cast<RealType>(36.0) * s2 +
                    static_cast<RealType>(16.0) * s4;
                if(isinf(f1) || isinf(cc))
                    return zero<RealType>;
                return f1*m_exp(-y*cc)/cc/three<RealType>;},
            zero<RealType>)/m_sqrt(three<RealType>);
    }

    /**
    * This function computes the integrand for the
    * G function (see validation script). It is not usable
    * on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_part the chi parameter of the particle
    * @param[in] csi the chi_photon/chi_part ratio
    *
    * @return the value of the integrand of the G function
    */
    template<typename RealType>
    RealType compute_G_integrand(
        const RealType chi_part, const RealType csi) noexcept
    {
        using namespace math;
        if( csi >= one<RealType> || chi_part == zero<RealType> )
            return zero<RealType>;

        if (csi == zero<RealType>)
            return std::numeric_limits<RealType>::infinity();

        const auto yy = compute_y(chi_part, csi);
        const RealType coeff = m_sqrt(three<>)/(two<>*pi<>);
        const auto inner = inner_integral(yy);

        const auto second_part = (csi*csi/(one<RealType>-csi))*
            k_v(two_thirds<RealType>,yy);
        if(std::isinf(second_part))
            return zero<RealType>;

        return coeff*csi*(inner + second_part);
    }

    /**
    * Computes the G function (see validation script).
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_part the chi parameter of the particle
    *
    * @return the value of the G function
    */
    template<typename RealType>
    RealType compute_G_function(const RealType chi_part)
    {
        using namespace math;
        return quad_a_b_s<RealType>(
            [=](RealType csi){
                return compute_G_integrand<RealType>(chi_part, csi)/csi;},
                zero<RealType>, one<RealType>);
    }

    /**
    * Auxiliary function to compute the numerator for the cumulative
    * probability distribution (see validation script).
    * It performs the integration from chi_start to chi_end.
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_particle the chi parameter of the particle
    * @param[in] chi_photon_start starting point of the integral
    * @param[in] chi_photon_end end point of the integral
    *
    * @return the value of the numerator of the cumulative probability distribution
    */
    template<typename RealType>
    RealType compute_cumulative_prob_numerator_a_b(
        const RealType chi_particle,
        RealType chi_photon_start,
        RealType chi_photon_end)
    {
        using namespace math;

        if(chi_particle  <= math::zero<RealType>) return math::zero<RealType>;
        if(chi_photon_end <= math::zero<RealType>) return math::zero<RealType>;

        if(chi_photon_end > chi_particle) chi_photon_end = chi_particle;
        if(chi_photon_start >= chi_photon_end) return zero<RealType>;

        auto frac_end = chi_photon_end/chi_particle;
        if(frac_end > math::one<RealType>) frac_end =  math::one<RealType>;
        auto frac_start = chi_photon_start/chi_particle;

        return math::quad_a_b_s<RealType>(
            [=](RealType csi){
                return compute_G_integrand<RealType>(chi_particle, csi)/csi;},
                frac_start, frac_end);
    }

    /**
    * Computes the numerator for the cumulative
    * probability distribution (see validation script).
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_particle the chi parameter of the particle
    * @param[in] chi_photon the chi parameter of the photon
    *
    * @return the value of the numerator of the cumulative probability distribution
    */
    template<typename RealType>
    RealType compute_cumulative_prob_numerator(
        const RealType chi_particle, RealType chi_photon)
    {
        using namespace math;
        return compute_cumulative_prob_numerator_a_b(
            chi_particle, zero<RealType>, chi_photon);
    }

    /**
    * Computes the cumulative probability distribution
    * (see validation script) for a vector of chi paramters
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_particle the chi parameter of the particle
    * @param[in] chi_photons the chi parameters of the photons
    *
    * @return the cumulative probability distribution calculated for all the chi parameters
    */
    template<typename RealType, typename VectorType>
    VectorType compute_cumulative_prob(
        const RealType chi_particle, const VectorType& chi_photons)
    {
        const auto den = compute_G_function(chi_particle);
        auto res = VectorType(chi_photons.size());

        if(chi_particle <= math::zero<RealType> || den <= math::zero<RealType>){
            for (auto& el: res ) el = math::zero<RealType>;
            return res;
        }

        std::transform(chi_photons.begin(), chi_photons.end(),
            res.begin(), [=](auto chi_phot){
                const auto val =
                    compute_cumulative_prob_numerator(chi_particle, chi_phot)/den;
                if(val <= math::one<RealType>) return val;
                return math::one<RealType>;});

        return res;
    }

    /**
    * Computes the cumulative probability distribution
    * (see validation script) for a vector of chi parameters.
    * Variant optimized for an ordered sequence of chi_photons.
    * It is not usable on GPUs.
    *
    * @tparam RealType the floating point type to be used
    *
    * @param[in] chi_phot the chi parameter of the particle
    * @param[in] chis the chi parameters of the photons (must be sorted)
    *
    * @return the cumulative probability distribution calculated for all the chi parameters
    */
    template<typename RealType, typename VectorType>
    VectorType compute_cumulative_prob_opt(
        const RealType chi_particle, const VectorType& chi_photons)
    {
        if(!std::is_sorted(chi_photons.begin(), chi_photons.end()))
            throw("Chi vector is not sorted!");

        using namespace math;
        const auto den = compute_G_function(chi_particle);
        auto res = VectorType(chi_photons.size());

        if(chi_particle <= zero<RealType> || den <= zero<RealType>){
            for (auto& el: res ) el = zero<RealType>;
            return res;
        }

        RealType old_chi = zero<RealType>;
        RealType sum = zero<RealType>;
        RealType c =  zero<RealType>;
        for (int i = 0; i < chi_photons.size(); ++i){
            //Kahan summation algorithm is used
            const auto y = compute_cumulative_prob_numerator_a_b(
                chi_particle, old_chi, chi_photons[i])/den - c;
            const auto t = sum + y;
            c = (t - sum) - y;
            sum = t;
            res[i] = sum;
            if(res[i] > one<RealType>) res[i] = one<RealType>;
            old_chi = chi_photons[i];
        }


        return res;
    }


}
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABULATED_FUNCTIONS
