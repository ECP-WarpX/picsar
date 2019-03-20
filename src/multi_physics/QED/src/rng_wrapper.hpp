#ifndef __PICSAR_MULTIPHYSICS_RNG_WRAPPER__
#define __PICSAR_MULTIPHYSICS_RNG_WRAPPER__

//This .hpp file provides a wrapper for the RNG.
//It always build a wrapper around the 64 bit Mersenne Twister availabe in the
// STL. If PXRMP_BUILD_WITH_KOKKOS_SUPPORT is defined, a support for Kokkos
//thread safe RNG is also built.

#include <random>
#include <cmath>
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

        //*** STL RNG Wrapper ***
        class stl_rng_wrapper
        {
        public:
            //The seed must be a 64 bit integer
            stl_rng_wrapper(uint64_t seed);
            //Alternatively, the constructor can take a generator (by move)
            stl_rng_wrapper(std::mt19937_64&& rng);

            //This line is just to make sure that copy and move constructor is generated
            stl_rng_wrapper(stl_rng_wrapper& ) = default;
            stl_rng_wrapper(stl_rng_wrapper&& ) = default;

            //Get rnd number uniformly distributed in [a,b)
            template<typename _REAL>
            PXRMP_FORCE_INLINE
            _REAL unf(_REAL a, _REAL b);

            //Get rnd number with exponential distribution
            template<typename _REAL>
            PXRMP_FORCE_INLINE
            _REAL exp(_REAL l);

        private:
            //Internally the 64 bit Mersenne Twister generator is used
            std::mt19937_64 rng;

        };

        //*** Kokkos RNG Wrapper ***
        //Only if built with the appropriate flag
        #ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

        //TO DO!

        //Kokkos-based wrapper is templated with respect to the generator pool
        //and the
        template<class generator_pool>
        class kokkos_rng_wrapper
        {
            public:
                //Constructor asks for a generator pool
                kokkos_rng_wrapper(generator_pool pool);

                //Copy constructor (this time explicitly defined to make Kokkos happy)
                kokkos_rng_wrapper(const kokkos_rng_wrapper& other);
                //Move constructor (this time explicitly defined to make Kokkos happy)
                kokkos_rng_wrapper(const kokkos_rng_wrapper&& other);

                //Get rnd number uniformly distributed in [a,b)
                template<typename _REAL>
                PXRMP_FORCE_INLINE
                _REAL unf(_REAL a, _REAL b) const;

                //Get rnd number with exponential distribution
                template<typename _REAL>
                PXRMP_FORCE_INLINE
                _REAL exp(_REAL l) const;

                //For Kokkos we need to
                //specialize for double and floats :-(
                //Non-template functions have precedence.
                PXRMP_FORCE_INLINE
                double unf(double a, double b) const;

                PXRMP_FORCE_INLINE
                float unf(float a, float b) const;

            private:
                //Internally, it keeps a copy of the generator pool
                const generator_pool pool;

                //Transforms a random number between [0.0, 1) into
                //a number extracted from exponential distribution
                template<typename _REAL>
                PXRMP_FORCE_INLINE
                _REAL unf2exp(_REAL uu) const;
        };

        #endif

    }
}

//############################################### Implementation

//*** STL RNG Wrapper ***

//Constructor with seed
picsar::multi_physics::stl_rng_wrapper::stl_rng_wrapper(uint64_t seed)
{
    rng.seed(seed);
}

//Constructor with move of an existing RNGt
picsar::multi_physics::stl_rng_wrapper::stl_rng_wrapper(std::mt19937_64&& rng):
    rng(rng){}

//Get rnd number uniformly distributed in [a,b)
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::stl_rng_wrapper::unf(_REAL a, _REAL b)
{
    auto unf_dist_a_b = std::uniform_real_distribution<_REAL>(a, b);
    return unf_dist_a_b(rng);
}

//Get rnd number with exponential distribution (l is the lambda param.)
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::stl_rng_wrapper::exp(_REAL l)
{
    auto exp_dist_l = std::exponential_distribution<_REAL>(l);
    return exp_dist_l(rng);
}


//*** Kokkos RNG Wrapper ***
//Only if built with the appropriate flag
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

//Constructor asks for a  pointer to a generator pool
template<class generator_pool>
picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
kokkos_rng_wrapper(generator_pool pool):
    pool{pool}
{}

//Copy constructor (this time explicitly defined to make Kokkos happy)
template<class generator_pool>
picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
kokkos_rng_wrapper(const kokkos_rng_wrapper& other):
    pool{other.pool}
{}

//Move constructor (this time explicitly defined to make Kokkos happy)
template<class generator_pool>
picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
kokkos_rng_wrapper(const kokkos_rng_wrapper&& other):
    pool{std::move(other.pool)}
{}

//Get rnd number uniformly distributed in [a,b)
//generic (to be specialized)
template<class generator_pool>
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
unf(_REAL a, _REAL b) const
{
    auto rand_gen = pool.get_state();

    //Since Kokkos RNG generates numbers between (a,b],
    //here we ask fot r in (-b,-a] and then we return -r
    _REAL res = -static_cast<_REAL>(rand_gen.drand(static_cast<_REAL>(-b),
        static_cast<_REAL>(-a)));
    pool.free_state(rand_gen);
    return res;
}

//Transforms a random number between [0.0, 1) into
//a number extracted from exponential distribution

template<class generator_pool>
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
 unf2exp(_REAL uu) const
{
    return -log(static_cast<_REAL>(1.0) - uu);
}

//Get rnd number with exponential distribution
//generic, should be specialized
template<class generator_pool>
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
exp(_REAL l) const {
    return unf2exp<_REAL>
        (unf(static_cast<_REAL>(0.0), static_cast<_REAL>(1.0)))/l;
}


//The specialized functions

template<class generator_pool>
PXRMP_FORCE_INLINE
double picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
unf(double a, double b) const
{
    auto rand_gen = pool.get_state();
    double res = -rand_gen.drand(-b, -a);
    pool.free_state(rand_gen);
    return res;
}

template<class generator_pool>
PXRMP_FORCE_INLINE
float picsar::multi_physics::kokkos_rng_wrapper<generator_pool>::
unf(float a, float b) const
{
    auto rand_gen = pool.get_state();
    float res = -rand_gen.frand(-b, -a);
    pool.free_state(rand_gen);
    return res;
}

#endif




#endif //__PICSAR_MULTIPHYSICS_RNG_WRAPPER__
