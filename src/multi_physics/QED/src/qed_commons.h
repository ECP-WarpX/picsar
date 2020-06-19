#ifndef PICSAR_MULTIPHYSICS_QED_COMMONS
#define PICSAR_MULTIPHYSICS_QED_COMMONS

/**
 * This header file should be included by all the source files of the library,
 * since it containes some common definitions.
 */


/**
 * The core functions of the library should run also on GPUs.
 * In order to do so, the user has to define PXRMP_GPU as follows:
 * #define PXRMP_GPU __host__ __device__
 * before including any file of the library.
 */
#ifndef PXRMP_GPU
  #define PXRMP_INTERNAL_NO_GPU
  #define PXRMP_INTERNAL_GPU_DECORATOR
#else
  #define PXRMP_INTERNAL_WITH_GPU
  #define PXRMP_INTERNAL_GPU_DECORATOR PXRMP_GPU
  #define PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY
#endif

/**
 * If GPU support is not enabled, there should be no need to use
 * GPU friendly PICSAR arrays. However, the user can force their
 * use with:
 * #define PXRMP_FORCE_PICSAR_ARRAY
 * This should be done only for debug purposes.
 */
#ifdef PXRMP_FORCE_PICSAR_ARRAY
   #define PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY
#endif

/**
 * The user can explicitly define a restrict keyword by doing
 * #define PXRMP_RESTRICT __restrict__
 * otherwise a choice based on the operating system is made.
 */
#ifdef PXRMP_RESTRICT
  #define PXRMR_INTERNAL_RESTRICT PXRMP_RESTRICT
#else
  #ifdef _WIN32
    #define PXRMR_INTERNAL_RESTRICT __restrict
  #else
    #define PXRMR_INTERNAL_RESTRICT __restrict__
  #endif
#endif

/**
 * The user can explicitly define a force inline keyword by doing
 * #define PXRMP_FORCE_INLINE __forceinline__
 * otherwise a choice based on the compiler is made
 */
#ifdef PXRMP_FORCE_INLINE
  #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR PXRMP_FORCE_INLINE
#else
  #if defined(__CUDA_ARCH__)
    #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR __forceinline__
  #elif defined(__INTEL_COMPILER)
    #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR inline __attribute__((always_inline))
  #elif defined(__clang__)
    #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR inline __attribute__((always_inline))
  #elif defined(__GNUC__)
    #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR inline __attribute__((always_inline))
  #elif defined(__ibmxl__)
    #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR inline __attribute__((always_inline))
  #else
    #define PXRMP_INTERNAL_FORCE_INLINE_DECORATOR inline
  #endif
#endif


/**
 * By doing
 * #define PXRMP_FORCE_BOOST_FOR_SPECFUNC
 * the user can force the use of Boost library for
 * special functions. Otherwise, a choice based on
 * the C++ language standard in use is made.
 * C++17 has several special functions defined in its standard
 * library.
 */
  #ifdef PXRMP_FORCE_BOOST_FOR_SPECFUNC
    #define PXRMP_INTERNAL_SPECFUNC_WITH_BOOST
  #else
    #if __cplusplus > 201402L
      #define PXRMP_INTERNAL_SPECFUNC_WITH_CXX17
    #else
      #define PXRMP_INTERNAL_SPECFUNC_WITH_BOOST
    #endif
  #endif

/**
 * By doing
 * #define PXRMP_FORCE_PICSAR_UPPER_BOUND
 * the user can force the use of a GPU-friendly PICSAR
 * implementation of the upper_bound algorithm.
 * Otherwise, a choice is made based on if GPU
 * support is enabled or not.
 */
  #ifdef PXRMP_FORCE_PICSAR_UPPER_BOUND
    #define PXRMP_INTERNAL_PICSAR_UPPER_BOUND
  #else
    #ifdef PXRMP_INTERNAL_WITH_GPU
      #define PXRMP_INTERNAL_PICSAR_UPPER_BOUND
    #else
      #define PXRMP_INTERNAL_STL_UPPER_BOUND
    #endif
  #endif

/**
 * If possible (i.e. if C++17 or more recent is used)
 * picsar makes use of "if constexpr". Otherwise, the
 * expression falls back to a regular "if".
 */
  #if __cplusplus > 201402L
    #define PXRMP_INTERNAL_CONSTEXPR_IF if constexpr
  #else
    #define PXRMP_INTERNAL_CONSTEXPR_IF if
  #endif

#endif// PICSAR_MULTIPHYSICS_QED_COMMONS
