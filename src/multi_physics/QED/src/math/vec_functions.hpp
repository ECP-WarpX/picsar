#ifndef PICSAR_MULTIPHYSICS_VEC_FUNCTIONS
#define PICSAR_MULTIPHYSICS_VEC_FUNCTIONS

#include <cmath>

//Should be included by all the src files of the library
#include "../qed_commons.h"

//Uses GPU-friendly arrays
#include "../containers/picsar_array.hpp"

namespace picsar{
namespace multi_physics{
namespace math{

    //This .hpp file contains functions to perform operations on 3-vectors
    // (norm, scalar multiplication, vector and cross product...)

    template <typename RealType>
    using vec3 = containers::picsar_array<RealType, 3>;

    /**
    * This function returns the squared norm of a 3-vector
    *
    * @tparam RealType the floating point type to be used
    * @param[in] the vector for which the squared norm has to be calculated
    * @return the squared norm
    */
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType norm2(const vec3<RealType>& vec) noexcept
    {
        return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
    }

    /**
    * This function returns the norm of a 3-vector
    *
    * @tparam RealType the floating point type to be used
    * @param[in] the vector for which the norm has to be calculated
    * @return the norm
    */
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType norm(const vec3<RealType>& vec) noexcept
    {
        return sqrt(norm2(vec));
    }

    /**
    * This function returns the dot product between two 3-vectors
    *
    * @tparam RealType the floating point type to be used
    * @param[in] va the first vector
    * @param[in] vb the second vector
    * @return the dot product between va and vb
    */
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dot(const vec3<RealType>& va, const vec3<RealType>& vb) noexcept
    {
        return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
    }

    /**
    * This function returns the cross product between two 3-vectors
    *
    * @tparam RealType the floating point type to be used
    * @param[in] va the first vector
    * @param[in] vb the second vector
    * @return the cross product between va and vb
    */
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

    /**
    * This function overloads the * operator for a vector and a scalar
    *
    * @tparam RealType the floating point type to be used
    * @param[in] v the vector
    * @param[in] s the scalar
    * @return the scalar product between v and s
    */
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

    /**
    * This function overloads the * operator for a scalar and a vector
    *
    * @tparam RealType the floating point type to be used
    * @param[in] s the scalar
    * @param[in] v the vector
    * @return the scalar product between s and v
    */
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

    /**
    * This function overloads the / operator for a vector and a scalar
    *
    * @tparam RealType the floating point type to be used
    * @param[in] v the vector
    * @param[in] s the scalar
    * @return v/s
    */
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

    /**
    * This function overloads the + operator for two vectors
    *
    * @tparam RealType the floating point type to be used
    * @param[in] va the first vector
    * @param[in] vb the second vector
    * @return va+vb
    */
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

    /**
    * This function overloads the - operator for two vectors
    *
    * @tparam RealType the floating point type to be used
    * @param[in] va the first vector
    * @param[in] vb the second vector
    * @return va-vb
    */
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
