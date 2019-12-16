#ifndef PICSAR_MULTIPHYSICS_POISSON_DISTRIB
#define PICSAR_MULTIPHYSICS_POISSON_DISTRIB

//size_t is defined here
#include <cstddef>

//Should be included by all the src files of the library
#include "../qed_commons.h"

namespace picsar{
namespace multi_physics{
namespace math{

//A function providing values extracted from a poisson distribution,
//given lambda and a number in the interval [0,1)
    template<typename T>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    size_t poisson_distrib
    (T lambda, T unf_zero_one_minus_epsi)
    {
        size_t k = 0;
        T p = exp(-lambda);
        T s = p;
        while (unf_zero_one_minus_epsi > s){
            const auto old_s = s;
            p = p*lambda/(++k);
            s += p;
            //If this is true we have reached the limit of the floating
            //point number that we are using
            if(s <= old_s)
                break;
        }
        return k;
    }

}
}
}


#endif //PICSAR_MULTIPHYSICS_POISSON_DISTRIB