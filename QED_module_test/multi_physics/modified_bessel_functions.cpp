#include "modified_bessel_functions.h"

using namespace picsar;

//C++17 has built-in support for modified bessel functions
#if __cplusplus > 201402L
double picsar::k_1_3(double){
    return std::cyl_bessel_k(1.0/3.0, x);
}

double picsar::k_2_3(double){
    return std::cyl_bessel_k(2.0/3.0, x);
}
#endif

double picsar::k_1_3(double){
    // TO IMPLEMENT
    return 0.0;
}

double picsar::k_2_3(double){
    // TO IMPLEMENT
    return 0.0;
}
