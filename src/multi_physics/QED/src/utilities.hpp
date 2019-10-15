#ifndef __PICSAR_MULTIPHYSICS_UTILITIES__
#define __PICSAR_MULTIPHYSICS_UTILITIES__

//This .hpp file contains general purpose functions to perform useful
//operations

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <limits>
#include <utility>

#ifndef PXRMP_CORE_ONLY
    //Uses the root finding algorithms provided by boost
    #include <boost/math/tools/roots.hpp>
#endif //PXRMP_CORE_ONLY

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //Generates a linearly spaced vector
        template<typename _REAL>
        std::vector<_REAL> generate_lin_spaced_vec
        (_REAL min, _REAL max, size_t size);

        //Generates a logarithmically spaced vector
        template<typename _REAL>
        std::vector<_REAL> generate_log_spaced_vec
        (_REAL min, _REAL max, size_t size);

        //Generates log_lin_log spaced vector
        template<typename _REAL>
        std::vector<_REAL> generate_log_lin_log_spaced_vec
        (_REAL min, _REAL max, size_t size);

        //GPU-friendly replacement of "std::upper_bound"
        template<typename T>
        PXRMP_GPU
        PXRMP_FORCE_INLINE
        const T* picsar_upper_bound(const T* first, const T* last, const T& val);

        //A function providing values extracted from a poisson distribution,
        //given lambda and a number in the interval [0,1)
        template<typename T>
        PXRMP_GPU
        PXRMP_FORCE_INLINE
        size_t poisson_distrib(T lambda, T unf_zero_one_minus_epsi);

#ifndef PXRMP_CORE_ONLY
        //A wrapper around the function provided by Boost library
        template<typename _REAL>
        _REAL bracket_and_solve_root
        (const std::function<_REAL(_REAL)>& f, _REAL guess, bool rising);
#endif
    }
}

//############################################### Implementation

template<typename _REAL>
std::vector<_REAL>
picsar::multi_physics::generate_log_spaced_vec
(_REAL min, _REAL max, size_t size)
{
    //Return empty vector upon error
    if( min >= max ||
        min < static_cast<_REAL>(0) || max< static_cast<_REAL>(0) ||
        size < 2)
        return std::vector<_REAL>(0);

    std::vector<_REAL> vec(size);

    _REAL val = min;
    _REAL mul = pow(max/min, static_cast<_REAL>(1.0/(size-1)));

    for(size_t i = 0; i < size-1; ++i){
        vec[i] = val;
        val*=mul;
    }
    vec.back() = max; //Enforces this exactly
    return vec;
}

//Generates a linearly spaced vector
template<typename _REAL>
std::vector<_REAL>
picsar::multi_physics::generate_lin_spaced_vec
(_REAL min, _REAL max, size_t size)
{
    //Return empty vector upon error
    if( min >= max || size < 2)
        return std::vector<_REAL>(0);

    std::vector<_REAL> vec(size);

    for(size_t i = 0; i < size-1; ++i){
        vec[i] = static_cast<_REAL>(i*(max-min)/(size-1.0)) + min;
    }
    vec.back() = max; //Enforces this exactly
    return vec;
}

//Generates log_lin_log spaced vector
template<typename _REAL>
std::vector<_REAL>
picsar::multi_physics::generate_log_lin_log_spaced_vec
(_REAL min, _REAL max, size_t size)
{
    std::vector<_REAL> vec(size);
    size_t size_first = size/3;
    size_t size_second = size/3;
    size_t size_third = size/3;

    _REAL first_val = max*static_cast<_REAL>(1.0/10.0);
    _REAL second_val = max*static_cast<_REAL>(9.0/10.0);

    size_t n = 0;
    std::generate(vec.begin(), vec.begin()+size_first,
    [=] () mutable { return min*exp((n++)*log(first_val/min)/(size_first)); });

    std::generate(vec.begin()+size_first, vec.begin()+size_first+size_second,
    [=] () mutable { return first_val + (second_val-first_val)*(n++)/(size_second); });

    std::generate(vec.begin()+size_first+size_second, vec.end(),
    [=] () mutable { return max*exp((size_third-1-(n++))*log(second_val/max)/(size_third-1)); });

    vec.front() = min;
    vec.back() = max;

    return vec;
}

template<typename T>
PXRMP_GPU
PXRMP_FORCE_INLINE
const T*
picsar::multi_physics::picsar_upper_bound
(const T* first, const T* last, const T& val)
{
    const T* it;
    size_t count, step;
    count = last-first;
    while(count>0){
        it = first;
        step = count/2;
        it += step;
         if (!(val<*it)){
             first = ++it;
             count -= step + 1;
         }
         else{
             count = step;
         }
    }
    return first;
}

//A function providing values extracted from a poisson distribution,
//given lambda and a number in the interval [0,1)
template<typename T>
PXRMP_GPU
PXRMP_FORCE_INLINE
size_t picsar::multi_physics::poisson_distrib
(T lambda, T unf_zero_one_minus_epsi)
{
    size_t k = 0;
    T p = exp(-lambda);
    T s = p;
    T old_s;
    while (unf_zero_one_minus_epsi > s){
        old_s = s;
        p = p*lambda/(++k);
        s += p;
        //If this is true we have reached the limit of the floating
        //point number that we are using
        if(s <= old_s)
            break;
    }
    return k;
}

#ifndef PXRMP_CORE_ONLY
    //A wrapper around the function provided by Boost library
    template<typename _REAL>
    _REAL picsar::multi_physics::bracket_and_solve_root
    (const std::function<_REAL(_REAL)>& f, _REAL guess, bool rising)
    {
        size_t digits = std::numeric_limits<_REAL>::digits;
        size_t precision_digits = digits - 2;
        boost::math::tools::eps_tolerance<_REAL> tol(precision_digits);

        _REAL factor = static_cast<_REAL>(2.0);

        size_t max_iter = 32;

        std::pair<_REAL, _REAL> r =
            boost::math::tools::bracket_and_solve_root
            (f, guess, factor, rising, tol, max_iter);

        return r.first + (r.second - r.first)/static_cast<_REAL>(2.0);
    }
#endif //PXRMP_CORE_ONLY

#endif // __PICSAR_MULTIPHYSICS_UTILITIES__
