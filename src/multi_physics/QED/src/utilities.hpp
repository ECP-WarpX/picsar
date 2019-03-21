#ifndef __PICSAR_MULTIPHYSICS_UTILITIES__
#define __PICSAR_MULTIPHYSICS_UTILITIES__

//This .hpp file contains general purpose functions to perform useful
//operations

#include <vector>
#include <cmath>
#include <functional>

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


#endif // __PICSAR_MULTIPHYSICS_UTILITIES__
