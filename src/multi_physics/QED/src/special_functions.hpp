#ifndef __PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS__
#define __PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS__

//This .hpp file is an extremely thin wrapper around special functions
//(Bessel functions for now) defined either in the STL (if C++17 is available)
//or in Boost library as a fallback.

//TODO: should we foresee a flag FORCE_USE_BOOST ?

//Set build option for the Bessel functions.
// 1) from STL (if C++ version > 14)
// 2) from Boost library
#if __cplusplus > 201402L
    #define PXRMP_SPECFUNC_WITH_CXX17
    #include <cmath>
#else
    #define PXRMP_SPECFUNC_WITH_BOOST
    #include <boost/math/special_functions/bessel.hpp>
#endif

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //For the moment we need just modified Bessel functions of the
        //second kind.
        //Different combinations of argument types can be accepted
        //(e.g. double + double or double + float).
    #ifdef PXRMP_SPECFUNC_WITH_CXX17
        template<typename _REAL_ARG1, typename _REAL_ARG2>
        auto k_v(_REAL_ARG1 v, _REAL_ARG2 x);
    #elif defined(PXRMP_SPECFUNC_WITH_BOOST)
        template<typename _REAL_ARG1, typename _REAL_ARG2>
        constexpr auto k_v(_REAL_ARG1 v, _REAL_ARG2 x)
        -> decltype(boost::math::cyl_bessel_k(v, x));
        //Ugly, ugly syntax. But unavoidable with C++11
    #endif

    }
}


//############################################### Implementation


#ifdef PXRMP_SPECFUNC_WITH_CXX17

template<typename _REAL_ARG1, typename _REAL_ARG2>
auto picsar::multi_physics::k_v(_REAL_ARG1 v, _REAL_ARG2 x)
{
    return std::cyl_bessel_k(v, x);
}

#elif defined(PXRMP_SPECFUNC_WITH_BOOST)

template<typename _REAL_ARG1, typename _REAL_ARG2>
constexpr auto picsar::multi_physics::k_v(_REAL_ARG1 v, _REAL_ARG2 x)
->decltype(boost::math::cyl_bessel_k(v, x))
{
    return boost::math::cyl_bessel_k(v, x);
}

#endif

#endif //__PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS__
