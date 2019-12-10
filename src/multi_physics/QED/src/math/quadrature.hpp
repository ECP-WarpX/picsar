#ifndef __PICSAR_MULTIPHYSICS_QUADRATURE__
#define __PICSAR_MULTIPHYSICS_QUADRATURE__

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

    //Provides integration in a finite interval [a,b]
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType quad_a_b(const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
        return boost::math::quadrature::trapezoidal(f, a, b);
    }

    //Provides integration in a finite interval [a,b] (more suitable for
    //singularities at boundaries )
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType quad_a_b_s(const std::function<RealType(RealType)>& f, RealType a, RealType b)
    {
        boost::math::quadrature::tanh_sinh<RealType> integrator;
        return integrator.integrate(f, a, b);
    }

    //Provides integration in a semi-infinite interval [a,inf)
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType quad_a_inf(const std::function<RealType(RealType)>& f, RealType a)
    {
        return boost::math::quadrature::gauss_kronrod<RealType, 15>
        ::integrate(f, a, std::numeric_limits<RealType>::infinity());
    }
}
}
}

#endif //__PICSAR_MULTIPHYSICS_QUADRATURE__
