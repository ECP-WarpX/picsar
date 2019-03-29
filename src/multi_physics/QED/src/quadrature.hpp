#ifndef __PICSAR_MULTIPHYSICS_QUADRATURE__
#define __PICSAR_MULTIPHYSICS_QUADRATURE__

//This .hpp file is a wrapper aroud the adaptive trapezoidal method
//provided by the Boost library
//NOTE: newer versions of Boost provide other (better) methods for
//quadrature

#include <cmath>
#include <functional>

#include <boost/math/quadrature/trapezoidal.hpp>

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
    return boost::math::quadrature::trapezoidal(f, a, b);
}

//This  function maps a semi-infinite interval to the finite interval [0,1],
//so that quad_a_b defined above could be used.
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quad_a_inf
(const std::function<_REAL(_REAL)>& f, _REAL a)
{
    _REAL zer = static_cast<_REAL>(0.0);
    _REAL one = static_cast<_REAL>(1.0);

    //Mapping function [0,1] --> [a, inf), defined as a lambda
    auto map_func = [one, a](_REAL x)
        { return a - one +  one/ (one - x);};

    //Analytical derivative of the mapping function
    auto deriv_map_func = [one](_REAL x)
        { return one/((one - x)*(one - x));};

    //Function to be integrated with  quad_a_b
    auto mod_func = [&f, &map_func, &deriv_map_func](_REAL x)
        { return f(map_func(x))*deriv_map_func(x);};

    //Alternative
    //
    // //Mapping function [0,1] --> [a, inf), defined as a lambda
    // auto map_func = [one, a](_REAL x)
    //     { return a - one - x + one/(one - x);};
    //
    // //Analytical derivative of the mapping function
    // auto deriv_map_func = [one](_REAL x)
    //     { return one/((one - x)*(one - x)) - one;};

    //WARNING: sqrt(numeric_limits<_REAL>::epsilon()) is added here to avoid
    //infinities. However, I am not very happy with this solution.
    //Suggestions?
    return picsar::multi_physics::quad_a_b<_REAL>
    (mod_func, zer, one-sqrt(std::numeric_limits<_REAL>::epsilon()));
}

#endif //__PICSAR_MULTIPHYSICS_QUADRATURE__
