#ifndef __PICSAR_MULTIPHYSICS_VEC_FUNCTIONS__
#define __PICSAR_MULTIPHYSICS_VEC_FUNCTIONS__

//This .hpp file contais functions to perform operations on 3-vectors
// (norm, scalar multiplication, vector and cross product...)

#include <cmath>
#include <array>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //A 3-vector
        template <typename _REAL>
        using vec3 = std::array<_REAL, 3>;

        //Norm of a 3-vector
        template <typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL norm(const vec3<_REAL>& vec);

        //Dot product of two 3-vectors
        template <typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL dot(const vec3<_REAL>& va, const vec3<_REAL>& vb);

        //Cross product of two 3-vectors
        template <typename _REAL>
        PXRMP_FORCE_INLINE
        vec3<_REAL> cross(const vec3<_REAL>& va, const vec3<_REAL>& vb);

    }
}

//############################################### Implementation

//Norm of a 3-vector
template <typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::norm(const vec3<_REAL>& vec)
{
    return sqrt(dot(vec, vec));
}

//Dot product of two 3-vectors
template <typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::dot(const vec3<_REAL>& va, const vec3<_REAL>& vb)
{
    return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
}

//Cross product of two 3-vectors
template <typename _REAL>
PXRMP_FORCE_INLINE
picsar::multi_physics::vec3<_REAL>
picsar::multi_physics::cross(const vec3<_REAL>& va, const vec3<_REAL>& vb)
{
    picsar::multi_physics::vec3<_REAL> out;
    out[0] = va[1]*vb[2] - va[2]*vb[1];
    out[1] = va[2]*vb[0] - va[0]*vb[2];
    out[2] = va[0]*vb[1] - va[1]*vb[0];
    return out;
}

#endif // __PICSAR_MULTIPHYSICS_VEC_FUNCTIONS__
