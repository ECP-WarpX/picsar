#ifndef PICSAR_MULTIPHYSICS_QUADRATURE
#define PICSAR_MULTIPHYSICS_QUADRATURE

#include <cmath>
#include <functional>
#include <limits>

#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

//Should be included by all the src files of the library
#include "../qed_commons.h"

namespace picsar{
namespace multi_physics{
namespace math{
    /**
    * This module is a wrapper is a wrapper aroud the trapezoidal,
    * gauss_kronrod, tanh_sinh & exp_sinh quadrature methods provided by the Boost library.
    * All the functions provided here accept a quadrature_algorithm template
    * parameter to choose which method will be used.
    */
    enum quadrature_algorithm {
        trapezoidal,
        tanh_sinh,
        exp_sinh,
        gauss_kronrod15,
        gauss_kronrod31,
        gauss_kronrod41,
        gauss_kronrod51,
        gauss_kronrod61
    };

    /**
    * This function performs the integration of the function f(x)
    * in the interval (a,b) using the method specified in the template parameter
    *
    * @tparam RealType the floating point type to be used
    * @tparam QuadAlgo the quadrature method to be used
    * @param[in] f the function which should be integrated
    * @param[in] a the left boundary of the integration region
    * @param[in] b the right boundary of the integration region
    * @return the integral of f in (a,b)
    */
    template<
        typename RealType, quadrature_algorithm QuadAlgo>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType generic_quad_a_b(
        const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
        if(QuadAlgo == quadrature_algorithm::trapezoidal){
            return boost::math::quadrature::trapezoidal(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::tanh_sinh){
            boost::math::quadrature::tanh_sinh<RealType> integrator;
            return integrator.integrate(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::exp_sinh){
            boost::math::quadrature::exp_sinh<RealType> integrator;
            return integrator.integrate(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::gauss_kronrod15){
            return boost::math::quadrature::gauss_kronrod<RealType, 15>
                ::integrate(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::gauss_kronrod31){
            return boost::math::quadrature::gauss_kronrod<RealType, 31>
                ::integrate(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::gauss_kronrod41){
            return boost::math::quadrature::gauss_kronrod<RealType, 41>
                ::integrate(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::gauss_kronrod51){
            return boost::math::quadrature::gauss_kronrod<RealType, 51>
                ::integrate(f, a, b);
        }
        else if(QuadAlgo == quadrature_algorithm::gauss_kronrod61){
            return boost::math::quadrature::gauss_kronrod<RealType, 61>
                ::integrate(f, a, b);
        }
        else
        {
            return boost::math::quadrature::trapezoidal(f, a, b);
        }
    }

    /**
    * This function performs the integration of the function f(x)
    * in the finite interval (a,b), using the "trapezoidal" quadrature method.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] f the function which should be integrated
    * @param[in] a the left boundary of the integration region
    * @param[in] b the right boundary of the integration region
    * @return the integral of f in (a,b)
    */
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType quad_a_b(
        const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
        return generic_quad_a_b<
            RealType, quadrature_algorithm::gauss_kronrod61>(f, a, b);
    }

    /**
    * This function performs the integration of the function f(x)
    * in the finite interval (a,b) using the "tanh_sinh" quadrature method,
    * to deal with possibile singularities at the boundaries.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] f the function which should be integrated
    * @param[in] a the left boundary of the integration region
    * @param[in] b the right boundary of the integration region
    * @return the integral of f in (a,b)
    */
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType quad_a_b_s(
        const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
         return generic_quad_a_b<RealType,
            quadrature_algorithm::tanh_sinh>(f, a, b);
    }

    /**
    * This function performs the integration of the function f(x)
    * in the semi-infinite interval (a,inf) using the "exp_sinh" quadrature method.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] a the left boundary of the integration region
    * @param[in] b the right boundary of the integration region
    * @return the integral of f in (a,b)
    */
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType quad_a_inf(
        const std::function<RealType(RealType)>& f, RealType a)
    {
        return generic_quad_a_b<RealType, quadrature_algorithm::exp_sinh>(
            f, a, std::numeric_limits<RealType>::infinity());
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUADRATURE
