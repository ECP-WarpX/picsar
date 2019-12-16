#ifndef PICSAR_MULTIPHYSICS_ALGO
#define PICSAR_MULTIPHYSICS_ALGO

//Should be included by all the src files of the library
#include "../qed_commons.h"

#ifdef PXRMP_INTERNAL_STL_UPPER_BOUND
    #include <algorithm>
#endif

namespace picsar{
namespace multi_physics{
namespace utils{

template<typename ForwardIt, typename T>
PXRMP_INTERNAL_GPU_DECORATOR
PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
ForwardIt
picsar_upper_bound
(ForwardIt first, ForwardIt last, const T& val)
{

#ifdef PXRMP_INTERNAL_PICSAR_UPPER_BOUND

    size_t count = last-first;
    do{
        auto it = first;
        const auto step = count/2;
        it += step;
         if (!(val<*it)){
             first = ++it;
             count -= step + 1;
         }
         else{
             count = step;
         }
    }while(count>0);

    return first;

#elif defined(PXRMP_INTERNAL_STL_UPPER_BOUND)
    return std::upper_bound(first, last, val);
#else
    #error Wrong preprocessor variable for picsar_upper_bound
#endif

}


}
}
}

#endif //PICSAR_MULTIPHYSICS_ALGO
