#ifndef PICSAR_MULTIPHYSICS_ALGO
#define PICSAR_MULTIPHYSICS_ALGO

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

#ifndef PXRMP_PICSAR_UPPER_BOUND
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
* If the preprocessor variable PXRMP_PICSAR_UPPER_BOUND
* is defined, a GPU-compatible version is used. Otherwise,
* the upper_bound function of the standard template library is used.
* The choice should be made automatically by the library depending on if
* the code is compiled for GPU. However, the PICSAR implementation can
* be forced defining PXRMP_PICSAR_UPPER_BOUND.
*
* @tparam ForwardIt the iterator type
* @tparam T the type of 'val'
* @param[in] first a ForwardIt pointing to the first element of the container
* @param[in] last a ForwardIt pointing to the end of the container (i.e. beyond the last element)
* @param[in] val the value to use to find the upper bound
* @return a ForwardIt to the upper bound
*/
template<typename ForwardIt, typename T>
PXRMP_GPU_QUALIFIER
PXRMP_FORCE_INLINE
ForwardIt
picsar_upper_bound
(ForwardIt first, ForwardIt last, const T& val)
{
#ifdef PXRMP_PICSAR_UPPER_BOUND

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

#else
    return std::upper_bound(first, last, val);
#endif
}

/**
* This function returns the value of the
* the first element in a sequence which compares greater than val.
* If no element in the range compares greater than val, the function returns last.
* In order to be completely general, the sequence is given by a functor
* f(int i) --> T res, which should respect the property "res2 >= res1 if i2 >= i1".
*
* @tparam T the type of 'val'
* @tparam Functor the type of Functor
* @param[in] first the first index (an int)
* @param[in] last the last index (an int)
* @param[in] val the value to use to find the upper bound
* @param[in] functor a functor f(int i) --> T res, respecting "res2 >= res1 if i2 >= i1"
* @return a ForwardIt to the upper bound
*/
template<typename T, typename Functor>
PXRMP_GPU_QUALIFIER
PXRMP_FORCE_INLINE
int
picsar_upper_bound_functor
(int first, const int last, const T& val, Functor&& functor)
{
    int count = last-first;
    do{
        auto i = first;
        const auto step = count/2;
        i += step;
         if (!(val<functor(i))){
             first = ++i;
             count -= step + 1;
         }
         else{
             count = step;
         }
    }while(count>0);

    return first;
}

/** \brief Performs a linear interpolation
 *
 * Performs a linear interpolation at x given the 2 points
 * (x0, f0) and (x1, f1)
 *
 * @param[in] x0
 * @param[in] x1
 * @param[in] f0
 * @param[in] f1
 * @param[in] x
 * @return the interpolation result
 */
template<typename RealType>
PXRMP_GPU_QUALIFIER
PXRMP_FORCE_INLINE
RealType linear_interp(
    RealType x0, RealType x1,
    RealType f0, RealType f1, RealType x)
{
    return ((x1-x)*f0 + (x-x0)*f1)/(x1-x0);
}

/** \brief Performs a bilinear interpolation
 *
 * Performs a bilinear interpolation at (x,y) given the 4 points
 * (x0, y0, f00), (x0, y1, f01), (x1, y0, f10), (x1, y1, f11).
 *
 * @param[in] x0
 * @param[in] x1
 * @param[in] y0
 * @param[in] y1
 * @param[in] f00
 * @param[in] f01
 * @param[in] f10
 * @param[in] f11
 * @param[in] f1
 * @param[in] x
 * @param[in] y
 * @return the interpolation result
 */
template<typename RealType>
PXRMP_GPU_QUALIFIER
PXRMP_FORCE_INLINE
RealType bilinear_interp(
    RealType x0, RealType x1,
    RealType y0, RealType y1,
    RealType f00, RealType f01,
    RealType f10, RealType f11,
    RealType x, RealType y)
{
    const auto fx0 = linear_interp(x0, x1, f00, f10, x);
    const auto fx1 = linear_interp(x0, x1, f01, f11, x);
    return linear_interp(y0, y1, fx0, fx1, y);
}

}
}
}

#endif //PICSAR_MULTIPHYSICS_ALGO
