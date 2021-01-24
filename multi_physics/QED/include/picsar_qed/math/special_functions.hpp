#ifndef PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS
#define PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//This .hpp file is an extremely thin wrapper around special functions
//(Bessel functions for now) defined either in the STL (if C++17 is available)
//or in Boost library as a fallback.
#ifdef PXRMP_USE_CXX17_FOR_SPECIAL_FUNCTIONS
    #include <cmath>
#else
    #include <boost/math/special_functions/bessel.hpp>
#endif

namespace picsar{
namespace multi_physics{
namespace math{

    /**
    * This function is a wrapper around the Bessel function
    * of the second kind defined either in the STL (if C++17 is available)
    * or in Boost library as a fallback (not usable on GPUs).
    *
    * @tparam RealType the floating point type to be used
    * @param[in] v order of the function
    * @param[in] x argument of the function
    * @return K_v(x)
    */
    template<typename RealType>
    inline RealType k_v(RealType v, RealType x)
    {
#ifdef PXRMP_USE_CXX17_FOR_SPECIAL_FUNCTIONS
        return std::cyl_bessel_k(v, x);
#else
        return boost::math::cyl_bessel_k(v, x);
#endif
    }

}
}
}

#endif //PICSAR_MULTIPHYSICS_SPECIAL_FUNCTIONS
