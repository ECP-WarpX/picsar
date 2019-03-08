#ifndef __PMP_QUADRATURE__
#define __PMP_QUADRATURE__

#include <stdexcept>
#include <string>
#include <functional>

//For the moment quadrature is based on trapezoidal quadrature from boost

#if HAS_BOOST_MATH
    #define QUADRATURE_BUILD_WITH_BOOST //If possible use integration from boost library
#else
    #define QUADRATURE_BUILD_WITH_BUILTIN //If not fall-back onto a builtin implementation
#endif


#ifdef QUADRATURE_BUILD_WITH_BOOST
  #include <boost/math/quadrature/trapezoidal.hpp>
#endif

namespace picsar{
    namespace multi_physics{
        double quad_a_b(const std::function<double(double)>& func, double a, double b);

        double quad_a_inf(const std::function<double(double)>& func, double a);
    }
}

#endif
