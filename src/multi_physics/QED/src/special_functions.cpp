#include "special_functions.h"

using namespace picsar::multi_physics;

//###############################################Bessel functions

double picsar::multi_physics::k_1_3(double x){
    return  k_v(1.0/3.0, x);
}

double picsar::multi_physics::k_2_3(double x){
    return k_v(2.0/3.0, x);
}

double picsar::multi_physics::k_5_3(double x){
    return k_v(5.0/3.0, x);
}

//C++17 has built-in support for modified bessel functions. Use this as first option
#ifdef BESSEL_BUILD_WITH_STL

double  picsar::multi_physics::k_v(double v, double x){
    return std::cyl_bessel_k(v, x);
}

//If not, link to boost libray
#elif defined(BESSEL_BUILD_WITH_BOOST)
double  picsar::multi_physics::k_v(double v, double x){
    return boost::math::cyl_bessel_k(v, x);
}

//If everything fails, use directly Boost source files (TO DO)
#else

double picsar::multi_physics::k_v(double v, double x){
    // TO IMPLEMENT (fails for safety reasons)
    throw std::logic_error(std::string("Function not yet implemented"));
    return 0.0;
}

#endif

//###############################################End bessel functions
