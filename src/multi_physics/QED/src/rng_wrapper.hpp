#ifndef __PICSAR_MULTIPHYSICS_RNG_WRAPPER__
#define __PICSAR_MULTIPHYSICS_RNG_WRAPPER__

//This .hpp file provides a wrapper for RNG.
//If PXRMP_BUILD_WITH_KOKKOS_SUPPORT is defined, a support for Kokkos
//thread safe RNG is also built.

#include <random>
#include <cstdint>

//Should be included by all the src files of the library
#include "qed_commons.h"

//Kokkos library is included only if the appropriate flag is set
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT
    #include <Kokkos_Random.hpp>
#endif


//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //This is an abstract class, which will be inherited by STL wrapper and
        //by the Kokkos wrapper
        template<typename _REAL>
        class rng_wrapper{
        public:
            //Get rnd number uniformly distributed in [a,b)
            PXRMP_FORCE_INLINE
            virtual _REAL unf(_REAL a, _REAL b) = 0;

            //Get rnd numbers with exponential distribution (l is the lambda param.)
            PXRMP_FORCE_INLINE
            virtual _REAL expl(_REAL l) = 0;
        };

        //Derived class to encapsulate STL rng (NOT thread-safe!)
        template<typename _REAL>
        class stl_rng_wrapper: public rng_wrapper<_REAL>
        {
        public:
            //Seed must be a 64 bit integer
            stl_rng_wrapper(int64_t seed);
            //Alternatively, constructor can take a generator
            stl_rng_wrapper(std::mt19937_64&& rng);

            //This line is just to make sure that the move constructor is generated
            stl_rng_wrapper(stl_rng_wrapper&& ) = default;

            //These functions are actually implemented by this class
            PXRMP_FORCE_INLINE
            _REAL unf(_REAL a, _REAL b);
            PXRMP_FORCE_INLINE
            _REAL expl(_REAL l);

        private:
            //Internally a 64 bit Mersenne Twister generator is used
            std::mt19937_64 rng;

        };

    }
}

//############################################### Implementation

//*** STL RNG Wrapper ***

//Constructor with seed
template<typename _REAL>
picsar::multi_physics::stl_rng_wrapper<_REAL>::stl_rng_wrapper(int64_t seed)
{
    rng.seed(seed);
}

//Constructor with move of an existing RNG
template<typename _REAL>
picsar::multi_physics::stl_rng_wrapper<_REAL>::stl_rng_wrapper
(std::mt19937_64&& rng):
    rng(rng){}

//Get rnd number uniformly distributed in [a,b)
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::stl_rng_wrapper<_REAL>::unf(_REAL a, _REAL b)
{
    auto unf_dist_a_b = std::uniform_real_distribution<_REAL>(a, b);
    return unf_dist_a_b(rng);
}

//Get rnd numbers with exponential distribution (l is the lambda param.)
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::stl_rng_wrapper<_REAL>::expl(_REAL l)
{
    auto exp_dist_l = std::exponential_distribution<_REAL>(l);
    return exp_dist_l(rng);
}


//*** Kokkos RNG Wrapper ***
//Only if built with the appropriate flag
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT
    //TODO
#endif




#endif //__PICSAR_MULTIPHYSICS_RNG_WRAPPER__
