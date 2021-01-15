#ifndef PICSAR_MULTIPHYSICS_MATH_CONSTANTS
#define PICSAR_MULTIPHYSICS_MATH_CONSTANTS

//Should be included by all the src files of the library
#include "picsar/src/multi_physics/QED/src/qed_commons.h"

namespace picsar{
namespace multi_physics{
namespace math{

    // Mathemtatical constants
    template<typename RealType = double>
    constexpr RealType pi = RealType(3.14159265358979323846264338327950288);
    //________________________

    // Useful templated constants
    template<typename RealType = double>
    constexpr RealType zero = RealType(0.0);

    template<typename RealType = double>
    constexpr RealType half = RealType(0.5);

    template<typename RealType = double>
    constexpr RealType one = RealType(1.0);

    template<typename RealType = double>
    constexpr RealType two = RealType(2.0);

    template<typename RealType = double>
    constexpr RealType three = RealType(3.0);

    template<typename RealType = double>
    constexpr RealType four = RealType(4.0);

    template<typename RealType = double>
    constexpr RealType one_third = RealType(1.0/3.0);

    template<typename RealType = double>
    constexpr RealType two_thirds = RealType(2.0/3.0);

    template<typename RealType = double>
    constexpr RealType five_thirds = RealType(5.0/3.0);
    //________________________
}
}
}

#endif //PICSAR_MULTIPHYSICS_MATH_CONSTANTS
