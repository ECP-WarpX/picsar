#include "modified_bessel_functions.h"

using namespace picsar::multi_physics;

//C++17 has built-in support for modified bessel functions. Use this as first option
#ifdef BESSEL_BUILD_WITH_STL
double k_1_3(double){
    return std::cyl_bessel_k(1.0/3.0, x);
}

double k_2_3(double){
    return std::cyl_bessel_k(2.0/3.0, x);
}

//If not, link to boost libray
#elif defined(BESSEL_BUILD_WITH_BOOST)
double k_1_3(double x){
     boost::math::cyl_bessel_k(1.0/3.0, x);
    return 0.0;
}

double k_2_3(double x){
     boost::math::cyl_bessel_k(2.0/3.0, x);
    return 0.0;
}

//If everything fails, use directly Boost source files (TO DO)
#else

double k_1_3(double x){
    // TO IMPLEMENT (fails for safety reasons)
    throw std::logic_error("Function not yet implemented");
    return 0.0;
}

double k_2_3(double x){
    // TO IMPLEMENT (fails for safety reasons)
    throw std::logic_error("Function not yet implemented");
    return 0.0;
}

#endif
