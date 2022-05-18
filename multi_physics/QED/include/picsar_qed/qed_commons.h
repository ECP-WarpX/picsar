#ifndef PICSAR_MULTIPHYSICS_QED_COMMONS
#define PICSAR_MULTIPHYSICS_QED_COMMONS

/**
 * This header file should be included by all the source files of the library,
 * since it containes some common definitions.
 */


/**
 * The core functions of the library run also on GPUs.
 * In order to do so, the user has to define PXRMP_GPU as follows:
 * #define PXRMP_WITH_GPU
 * Otherwise the code is compiled without GPU support
 */
//#define PXRMP_WITH_GPU


/**
* If the code is compiled with GPU support,
 * PICSAR sets PXRMP_GPU_QUALIFIER to `__host__ __device__`
 * However, the user can override this by defining PXRMP_GPU_QUALIFIER
 * manually.
 */
#ifdef PXRMP_WITH_GPU
    #ifndef PXRMP_GPU_QUALIFIER
        #define PXRMP_GPU_QUALIFIER __host__ __device__
    #endif
#else
    #define PXRMP_GPU_QUALIFIER
#endif


/**
 * If GPU support is not enabled, there should be no need to use
 * GPU friendly PICSAR arrays. However, the user can force their
 * use with:
 * #define PXRMP_FORCE_PICSAR_ARRAY
 * This should be done only for debug purposes.
 */
#if defined(PXRMP_WITH_GPU) || defined(PXRMP_FORCE_PICSAR_ARRAY)
    #define PXRMP_ENABLE_GPU_FRIENDLY_ARRAY
#endif


/**
 * The user can explicitly define a restrict keyword by doing
 * #define PXRMP_RESTRICT __restrict__
 * otherwise a choice based on the operating system is made.
 */
#ifndef PXRMP_RESTRICT
  #ifdef _WIN32
    #define PXRMR_RESTRICT __restrict
  #else
    #define PXRMR_RESTRICT __restrict__
  #endif
#endif


/**
 * The user can explicitly define a force inline keyword by doing
 * #define PXRMP_FORCE_INLINE __forceinline__
 * otherwise a choice based on the compiler is made
 */
#ifndef PXRMP_FORCE_INLINE
  #if defined(__CUDA_ARCH__)
    #define PXRMP_FORCE_INLINE __forceinline__
  #elif defined(__INTEL_COMPILER)
    #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
  #elif defined(__clang__)
    #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
  #elif defined(__GNUC__)
    #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
  #elif defined(__ibmxl__)
    #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
  #else
    #define PXRMP_FORCE_INLINE inline
  #endif
#endif


/**
 * By doing
 * #define PXRMP_USE_CXX17_FOR_SPECIAL_FUNCTIONS
 * the user can choose to use special functions provided by the C++17
 * standard library if C++17 is available. Otherwise, special functions
 * provided by Boost are used.
 */
 #ifdef PXRMP_USE_CXX17_FOR_SPECIAL_FUNCTIONS
     #if __cplusplus < 201703L
        #error C++17 or above is needed to enable special functions from the standard C++ library
     #endif
 #endif


/**
 * If GPU support is enabled, the library uses an internal, GPU-friendly,
 * implementation of the upper_bound and of the lower_bound algorithm.
 * For test purposes, this implementation
 * can be used also on CPU, instead of the implementation provided by the STL.
 * This can be achieved by doing:
 * #define PXRMP_PICSAR_UPPER_BOUND
 * #define PXRMP_PICSAR_LOWER_BOUND
 */
#ifdef PXRMP_WITH_GPU
    #define PXRMP_PICSAR_UPPER_BOUND
    #define PXRMP_PICSAR_LOWER_BOUND
#endif


/**
 * If possible (i.e. if C++17 or more recent is used)
 * picsar makes use of "if constexpr". Otherwise, the
 * expression falls back to a regular "if".
 */
  #if __cplusplus >= 201703L
    #define PXRMP_CONSTEXPR_IF if constexpr
  #else
    #define PXRMP_CONSTEXPR_IF if
  #endif

 /**
 * Unless PXRMP_PREVENT_USE_STD_FOR_MATH is defined by the
 * user std::sqrt, std::cbrt... mathematical functions
 * are used.
 */
//#define PXRMP_PREVENT_USE_STD_FOR_MATH


/**
* PXRMP_DPCPP_FIX enables a workaround to allow compilation with
* DPC++, which apparently has issues with floorf or std::floor(x) when
* x is a float
*/
//
#ifdef __SYCL_DEVICE_ONLY__
#   ifndef PXRMP_DPCPP_FIX
#       define PXRMP_DPCPP_FIX 1
#   endif
#endif

/**
* PXRMP_HAS_OPENMP enables the use of OpenMP to calculate lookup tables
* on the CPU.
*/
//#define PXRMP_HAS_OPENMP

#endif// PICSAR_MULTIPHYSICS_QED_COMMONS
