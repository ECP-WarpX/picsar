#ifndef PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS
#define PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS

//This .hpp file is an extremely thin wrapper around special functions
//(Bessel functions for now) defined either in the STL (if C++17 is available)
//or in Boost library as a fallback.

//Should be included by all the src files of the library
#include "../qed_commons.h"

#ifdef PXRMP_INTERNAL_SPECFUNC_WITH_CXX17
    #include <cmath>
#elif defined(PXRMP_INTERNAL_SPECFUNC_WITH_BOOST)
    #include <boost/math/special_functions/bessel.hpp>
#endif

namespace picsar{
namespace multi_physics{
namespace math{

    /**
    * This function is a wrapper around the Bessel function
    * of the second kind defined either in the STL (if C++17 is available)
    * or in Boost library as a fallback.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] v order of the function
    * @param[in] x argument of the function
    * @return K_v(x)
    */
    template<typename RealType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType k_v(RealType v, RealType x)
    {
#ifdef PXRMP_INTERNAL_SPECFUNC_WITH_CXX17
        return std::cyl_bessel_k(v, x);
#elif defined(PXRMP_INTERNAL_SPECFUNC_WITH_BOOST)
        return boost::math::cyl_bessel_k(v, x);
#endif
    }

}
}
}

#endif //PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS
