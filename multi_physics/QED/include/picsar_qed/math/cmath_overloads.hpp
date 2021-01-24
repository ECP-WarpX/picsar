#ifndef PICSAR_MULTIPHYSICS_CMATH_OVERLOADS
#define PICSAR_MULTIPHYSICS_CMATH_OVERLOADS

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

#include <cmath>

#ifdef PXRMP_DPCPP_FIX
    #include <CL/sycl.hpp>
#endif

namespace picsar{
namespace multi_physics{
namespace math{

#ifdef PXRMP_PREVENT_USE_STD_FOR_MATH


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
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_sqrt(const RealType x) noexcept
    {
        return sqrt(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_sqrt(const double x) noexcept
    {
        return sqrt(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
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
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_cbrt(const RealType x) noexcept
    {
        return cbrt(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_cbrt(const double x) noexcept
    {
        return cbrt(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
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
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_log(const RealType x) noexcept
    {
        return log(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_log(const double x) noexcept
    {
        return log(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
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
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_exp(const RealType x) noexcept
    {
        return exp(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_exp(const double x) noexcept
    {
        return exp(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
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
    * @return tanh of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_tanh(const RealType x) noexcept
    {
        return tanh(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_tanh(const double x) noexcept
    {
        return tanh(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    float m_tanh(const float x) noexcept
    {
        return tanhf(x);
    }

    //________________________________________________________

    /**
    * This function replaces the overload of floor provided
    * by the Standard Template Library. It calls either
    * floor or floorf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the floor of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_floor(const RealType x) noexcept
    {
        return floor(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_floor(const double x) noexcept
    {
        return floor(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    float m_floor(const float x) noexcept
    {
#ifdef PXRMP_DPCPP_FIX
        return cl::sycl::floorf(x);
#else
        return floorf(x);
#endif
    }

    /**
    * This function replaces the overload of fabs provided
    * by the Standard Template Library. It calls either
    * fabs or fabsf in cmath depending on the variable type.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the fabs of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_fabs(const RealType x) noexcept
    {
        return fabs(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    double m_fabs(const double x) noexcept
    {
        return fabs(x);
    }

    template<>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    float m_fabs(const float x) noexcept
    {
        return fabsf(x);
    }

#else

    /**
    * This function wraps the overload of sqrt provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the square root of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_sqrt(const RealType x) noexcept
    {
        return std::sqrt(x);
    }

    /**
    * This function wraps the overload of cbrt provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the cubic root of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_cbrt(const RealType x) noexcept
    {
        return std::cbrt(x);
    }

    /**
    * This function wraps the overload of log provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the log of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_log(const RealType x) noexcept
    {
        return std::log(x);
    }

    /**
    * This function wraps the overload of exp provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the exponential of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_exp(const RealType x) noexcept
    {
        return std::exp(x);
    }

    /**
    * This function wraps the overload of tanh provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return tanh of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_tanh(const RealType x) noexcept
    {
        return std::tanh(x);
    }

    /**
    * This function wraps the overload of floor provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the floor of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_floor(const RealType x) noexcept
    {
#ifdef PXRMP_DPCPP_FIX
        return cl::sycl::floor(x);
#else
        return std::floor(x);
#endif
    }

    /**
    * This function wraps the overload of fabs provided
    * by the Standard Template Library.
    *
    * @tparam RealType the floating point type to be used
    * @param[in] x argument
    * @return the fabs of x
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
    RealType m_fabs(const RealType x) noexcept
    {
        return std::fabs(x);
    }

#endif

}
}
}


#endif //PICSAR_MULTIPHYSICS_CMATH_OVERLOADS
