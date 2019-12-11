#ifndef PICSAR_MULTIPHYSICS_VEC_FUNCTIONS
#define PICSAR_MULTIPHYSICS_VEC_FUNCTIONS

//This .hpp file contains functions to perform operations on 3-vectors
// (norm, scalar multiplication, vector and cross product...)

#include <cmath>

//Should be included by all the src files of the library
#include "../qed_commons.h"

//Uses GPU-friendly arrays
#include "../containers/picsar_array.hpp"

//############################################### Declaration

namespace picsar{
namespace multi_physics{
namespace math{

    //A 3-vector
    template <typename RealType>
    using vec3 = containers::picsar_array<RealType, 3>;

    //Squared norm of a 3-vector
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType norm2(const vec3<RealType>& vec) noexcept
    {
        return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
    }

    //Norm of a 3-vector
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType norm(const vec3<RealType>& vec) noexcept
    {
        return sqrt(norm2(vec));
    }

    //Dot product of two 3-vectors
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dot(const vec3<RealType>& va, const vec3<RealType>& vb) noexcept
    {
        return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
    }

    //Cross product of two 3-vectors
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    vec3<RealType> cross(const vec3<RealType>& va, const vec3<RealType>& vb) noexcept
    {
        vec3<RealType> out;
        out[0] = va[1]*vb[2] - va[2]*vb[1];
        out[1] = va[2]*vb[0] - va[0]*vb[2];
        out[2] = va[0]*vb[1] - va[1]*vb[0];
        return out;
    }

    //Product of a vector times a scalar
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    vec3<RealType> operator*(const vec3<RealType>& v, const RealType s) noexcept
    {
        vec3<RealType> out;
        out[0] = v[0]*s;
        out[1] = v[1]*s;
        out[2] = v[2]*s;
        return out;
    }

    //Product of a scalar times a vector
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    vec3<RealType> operator*(const RealType s, const vec3<RealType>& v) noexcept
    {
        vec3<RealType> out;
        out[0] = v[0]*s;
        out[1] = v[1]*s;
        out[2] = v[2]*s;
        return out;
    }

    //Vector divided by a scalar
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    vec3<RealType> operator/(const vec3<RealType>& v, const RealType s) noexcept
    {
        vec3<RealType> out;
        out[0] = v[0]/s;
        out[1] = v[1]/s;
        out[2] = v[2]/s;
        return out;
    }

    //Add two vectors
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    vec3<RealType> operator+(const vec3<RealType>& va, const vec3<RealType>& vb) noexcept
    {
        vec3<RealType> out;
        out[0] = va[0] + vb[0];
        out[1] = va[1] + vb[1];
        out[2] = va[2] + vb[2];
        return out;
    }

    //Subtract two vectors
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    vec3<RealType> operator-(const vec3<RealType>& va, const vec3<RealType>& vb) noexcept
    {
        vec3<RealType> out;
        out[0] = va[0] - vb[0];
        out[1] = va[1] - vb[1];
        out[2] = va[2] - vb[2];
        return out;
    }

}
}
}

#endif // PICSAR_MULTIPHYSICS_VEC_FUNCTIONS
