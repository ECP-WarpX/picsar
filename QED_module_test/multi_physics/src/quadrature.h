#ifndef __PMP_QUADRATURE__
#define __PMP_QUADRATURE__

#if HAS_BOOST_MATH
    #define QUADRATURE_BUILD_WITH_BOOST //If possible use integration from boost library
#else
    #define QUADRATURE_BUILD_WITH_BUILTIN //If not fall-back onto a builtin implementation
#endif


#ifdef QUADRATURE_BUILD_WITH_BOOST
//TO DO
#endif
//TO DO
#endif
