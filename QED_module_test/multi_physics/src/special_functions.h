#ifndef __PMP_SPEC_FUNCTIONS__
#define __PMP_SPEC_FUNCTIONS__

#include <cmath>

//Set build option for the Bessel functions.
// 1) from STL (if C++17)
// 2) from Boost library (if available)
// 3) using built-in code

#if __cplusplus > 201402L
    #define BESSEL_BUILD_WITH_STL
#else
    #if HAS_BOOST_MATH
        #define BESSEL_BUILD_WITH_BOOST
    #else
        #define BESSEL_BUILD_WITH_BUILTIN
    #endif
#endif


#ifdef BESSEL_BUILD_WITH_BOOST
    #include <boost/math/special_functions/bessel.hpp>
#endif

namespace picsar{
    namespace multi_physics{
        double k_1_3(double);

        double k_2_3(double);
    }
}

#endif
