#ifndef __PICSAR_MULTIPHYSICS_QUADRATURE__
#define __PICSAR_MULTIPHYSICS_QUADRATURE__

//This .hpp file is a wrapper aroud the tanh_sinh  method
//provided by the Boost library

#include <cmath>
#include <functional>
#include <limits>

#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //Provides integration in a finite interval [a,b]
        template<typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL quad_a_b(const std::function<_REAL(_REAL)>& f, _REAL a, _REAL b);

        //Provides integration in a semi-infinite interval [a,inf)
        template<typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL quad_a_inf(const std::function<_REAL(_REAL)>& f, _REAL a);
    }
}

//############################################### Implementation

//This is really just a wrapper around the function provided by
//the Boost libray
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quad_a_b
(const std::function<_REAL(_REAL)>& f, _REAL a, _REAL b)
{
    //    boost::math::quadrature::tanh_sinh<_REAL> integrator;
    //    return integrator.integrate(f, a, b);
    return boost::math::quadrature::trapezoidal(f, a, b);
}

//This is really just a wrapper around the function provided by
//the Boost libray
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quad_a_inf
(const std::function<_REAL(_REAL)>& f, _REAL a)
{
    boost::math::quadrature::exp_sinh<_REAL> integrator;
    return integrator.integrate(f, a, std::numeric_limits<_REAL>::infinity());
}

#endif //__PICSAR_MULTIPHYSICS_QUADRATURE__
