#ifndef PICSAR_MULTIPHYSICS_QUADRATURE
#define PICSAR_MULTIPHYSICS_QUADRATURE

//This .hpp file is a wrapper aroud the trapezoidal, gauss_kronrod & tanh_sinh methods
//provided by the Boost library

#include <cmath>
#include <functional>
#include <limits>

#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

//Should be included by all the src files of the library
#include "../qed_commons.h"

//############################################### Declaration

namespace picsar{
namespace multi_physics{
namespace math{

    enum quadrature_algorithm {
        trapezoidal,
        tanh_sinh,
        gauss_kronrod15,
        gauss_kronrod31,
        gauss_kronrod41,
        gauss_kronrod51,
        gauss_kronrod61
    };

    //Generic
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

    //Provides integration in a finite interval [a,b]
    template<
        typename RealType,
        quadrature_algorithm QuadAlgo = quadrature_algorithm::trapezoidal>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType quad_a_b(
        const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
        return generic_quad_a_b<RealType, QuadAlgo>(f, a, b);
    }

    //Provides integration in a finite interval [a,b] (more suitable for
    //singularities at boundaries )
    template<
        typename RealType,
        quadrature_algorithm QuadAlgo = quadrature_algorithm::tanh_sinh>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType quad_a_b_s(
        const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
         return generic_quad_a_b<RealType, QuadAlgo>(f, a, b);
    }

    //Provides integration in a semi-infinite interval [a,inf)
    template<
        typename RealType,
        quadrature_algorithm QuadAlgo = quadrature_algorithm::gauss_kronrod15>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType quad_a_inf(
        const std::function<RealType(RealType)>& f, RealType a)
    {
        return generic_quad_a_b<RealType, QuadAlgo>(
            f, a, std::numeric_limits<RealType>::infinity());
    }
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUADRATURE
