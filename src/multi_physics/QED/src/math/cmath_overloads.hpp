#ifndef PICSAR_MULTIPHYSICS_CMATH_OVERLOADS
#define PICSAR_MULTIPHYSICS_CMATH_OVERLOADS

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <cmath>

namespace picsar{
namespace multi_physics{
namespace math{

    /**
    * This function replaces the overload of sqrt provided
    * by the Standard Template Library. It calls either
    * sqrt or sqrtf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the square root of x
    */
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType m_sqrt(const RealType x) noexcept
    {
        return sqrt(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    double m_sqrt(const double x) noexcept
    {
        return sqrt(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    float m_sqrt(const float x) noexcept
    {
        return sqrtf(x);
    }

    //________________________________________________________

    /**
    * This function replaces the overload of cbrt provided
    * by the Standard Template Library. It calls either
    * cbrt or cbrtf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the cubic root of x
    */
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType m_cbrt(const RealType x) noexcept
    {
        return cbrt(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    double m_cbrt(const double x) noexcept
    {
        return cbrt(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    float m_cbrt(const float x) noexcept
    {
        return cbrtf(x);
    }

    //________________________________________________________

    /**
    * This function replaces the overload of log provided
    * by the Standard Template Library. It calls either
    * log or logf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the log of x
    */
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType m_log(const RealType x) noexcept
    {
        return log(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    double m_log(const double x) noexcept
    {
        return log(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    float m_log(const float x) noexcept
    {
        return logf(x);
    }

    //________________________________________________________

    /**
    * This function replaces the overload of exp provided
    * by the Standard Template Library. It calls either
    * exp or expf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the exponential of x
    */
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType m_exp(const RealType x) noexcept
    {
        return exp(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    double m_exp(const double x) noexcept
    {
        return exp(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    float m_exp(const float x) noexcept
    {
        return expf(x);
    }

    //________________________________________________________

    /**
    * This function replaces the overload of tanh provided
    * by the Standard Template Library. It calls either
    * tanh or tanhf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the exponential of x
    */
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType m_tanh(const RealType x) noexcept
    {
        return tanh(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    double m_tanh(const double x) noexcept
    {
        return tanh(x);
    }

    template<>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    float m_tanh(const float x) noexcept
    {
        return tanhf(x);
    }

}
}
}

#endif //PICSAR_MULTIPHYSICS_CMATH_OVERLOADS