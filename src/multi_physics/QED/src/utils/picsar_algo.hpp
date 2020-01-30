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

/**
* This function returns an iterator pointing
* to the first element in the range [first,last)
* which compares greater than val. If no element in the range compares
* greater than val, the function returns last.
* If the internal preprocessor variable PXRMP_INTERNAL_PICSAR_UPPER_BOUND
* is defined, a GPU-compatible version is used. Otherwise, if the
* internal preprocessor variable PXRMP_INTERNAL_STL_UPPER_BOUND is defined,
* the upper_bound function of the standard template library is used.
* The choice should be made automatically by the library depending on if
* the code is compiled for GPU. However, the PICSAR implementation can
* be forced defining PXRMP_FORCE_PICSAR_UPPER_BOUND.
*
* @tparam ForwardIt the iterator type
* @tparam the type of 'val'
* @param[in] first a ForwardIt pointing to the first element of the container
* @param[in] last a ForwardIt pointing to the end of the container (i.e. beyond the last element)
* @param[in] val the value to use to find the upper bound
* @return a ForwardIt to the upper bound
*/
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
