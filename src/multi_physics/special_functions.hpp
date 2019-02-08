#ifndef __SPECIAL_FUNCTIONS__
#define __SPECIAL_FUNCTIONS__

#include <cmath>

// Modified Bessel from library if C++17. If not, use implementation from
// numerical recipes.

#if __cplusplus==201703L
const auto& modified_bessel_K = static_cast<double>(cyl_bessel_k);
#else
double modified_bessel_K(double v, double x){
    //TODO
    return 0.0;
}
#endif

#endif
