#ifndef PICSAR_MULTIPHYSICS_RNG_WRAPPER
#define PICSAR_MULTIPHYSICS_RNG_WRAPPER

//This .hpp file provides a wrapper for the RNG.
//It builds a wrapper around the 64 bit Mersenne Twister availabe in the STL.

#include <random>
#include <cmath>
#include <cstdint>

//Should be included by all the src files of the library
#include "../qed_commons.h"

namespace picsar{
namespace multi_physics{
namespace utils{

    //*** STL RNG Wrapper ***
    class stl_rng_wrapper
    {
    public:
        //The seed must be a 64 bit integer
        stl_rng_wrapper(uint64_t seed):
            m_rng{seed}{}

        //Get rnd number uniformly distributed in [a,b)
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType unf(RealType a, RealType b)
        {
            auto unf_dist_a_b =
                std::uniform_real_distribution<RealType>(a, b);
            return unf_dist_a_b(m_rng);
        }

        //Get rnd number with exponential distribution
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType exp(RealType l)
        {
            auto exp_dist_l =
                std::exponential_distribution<RealType>(l);
            return exp_dist_l(m_rng);
        }

        //Get rnd number with poisson distribution
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        size_t poisson(RealType l)
        {
            auto poisson_dist_l =
                std::poisson_distribution<size_t>(l);
            return poisson_dist_l(m_rng);
        }

        //Get rnd number with gaussian distribution
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType gaussian(RealType mean, RealType deviation)
        {
            auto gaussian_dist_l =
                std::normal_distribution<RealType>(mean, deviation);
            return gaussian_dist_l(m_rng);
        }

    private:
        //Internally the 64 bit Mersenne Twister generator is used
        std::mt19937_64 m_rng;

    };

}
}
}

#endif //PICSAR_MULTIPHYSICS_RNG_WRAPPER
