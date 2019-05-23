#ifndef __PICSAR_MULTIPHYSICS_RNG_WRAPPER__
#define __PICSAR_MULTIPHYSICS_RNG_WRAPPER__

//This .hpp file provides a wrapper for the RNG.
//It builds a wrapper around the 64 bit Mersenne Twister availabe in the STL.

#include <random>
#include <cmath>
#include <cstdint>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration
namespace picsar{
    namespace multi_physics{

        //*** STL RNG Wrapper ***
        template<typename _REAL>
        class stl_rng_wrapper
        {
        public:
            //The seed must be a 64 bit integer
            stl_rng_wrapper(uint64_t seed);
            //Alternatively, the constructor can take a generator (by move)
            stl_rng_wrapper(std::mt19937_64&& rng);

            //This line is just to make sure that copy & move constructors
            //are generated
            stl_rng_wrapper(const stl_rng_wrapper& ) = default;
            stl_rng_wrapper(stl_rng_wrapper&& ) = default;

            //Assignment operator
            stl_rng_wrapper& operator= (const stl_rng_wrapper& other);

            //Get rnd number uniformly distributed in [a,b)
            PXRMP_FORCE_INLINE
            _REAL unf(_REAL a, _REAL b);

            //Get rnd number with exponential distribution
            PXRMP_FORCE_INLINE
            _REAL exp(_REAL l);

        private:
            //Internally the 64 bit Mersenne Twister generator is used
            std::mt19937_64 rng;

        };
    }
}

//############################################### Implementation

//*** STL RNG Wrapper ***

//Constructor with seed
template<typename _REAL>
picsar::multi_physics::stl_rng_wrapper<_REAL>::stl_rng_wrapper(uint64_t seed)
{
    rng.seed(seed);
}

 //Constructor with move of an existing RNG
 template<typename _REAL>
 picsar::multi_physics::stl_rng_wrapper<_REAL>::stl_rng_wrapper
 (std::mt19937_64&& rng):
    rng(std::move(rng)){}

//Assignment operator
template<typename _REAL>
picsar::multi_physics::stl_rng_wrapper<_REAL>&
picsar::multi_physics::stl_rng_wrapper<_REAL>::operator=
(const stl_rng_wrapper<_REAL>& other)
{
    rng = other.rng;
    return *this;
}

//Get rnd number uniformly distributed in [a,b)
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::stl_rng_wrapper<_REAL>::unf(_REAL a, _REAL b)
{
    auto unf_dist_a_b = std::uniform_real_distribution<_REAL>(a, b);
    return unf_dist_a_b(rng);
}

//Get rnd number with exponential distribution (l is the lambda param.)
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::stl_rng_wrapper<_REAL>::exp(_REAL l)
{
    auto exp_dist_l = std::exponential_distribution<_REAL>(l);
    return exp_dist_l(rng);
}

#endif //__PICSAR_MULTIPHYSICS_RNG_WRAPPER__
