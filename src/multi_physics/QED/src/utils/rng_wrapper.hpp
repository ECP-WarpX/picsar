#ifndef PICSAR_MULTIPHYSICS_RNG_WRAPPER
#define PICSAR_MULTIPHYSICS_RNG_WRAPPER

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <random>
#include <cmath>
#include <cstdint>

namespace picsar{
namespace multi_physics{
namespace utils{

    /**
    * This class is a wrapper around the 64 bit Mersenne Twister generator
    * and the distribution functions provided by the standard template library.
    * It can be used to instantiate the engine objects of the high order interface
    * of the library. However, any other used-defined class implementing the
    * same interface can be used as well.
    */
    class stl_rng_wrapper
    {
    public:

        /**
        * The constructor initializes the random number generator
        *
        * @param[in] seed a 64bit integer
        */
        stl_rng_wrapper(uint64_t seed):
            m_rng{seed}{}


        /**
        * This function returns floating point numbers uniformly
        * distributed in [a,b).
        *
        * @tparam RealType the floating point number type
        * @param[in] a the lower limit of the interval
        * @param[in] b the upper limit of the interval
        * @return a random number uniformly distributed in [a,b)
        */
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType unf(RealType a, RealType b)
        {
            auto unf_dist_a_b =
                std::uniform_real_distribution<RealType>(a, b);
            return unf_dist_a_b(m_rng);
        }

        /**
        * This function returns positive floating point numbers with
        * an exponential probability distribution f(x) = l * exp(-l * x)
        *
        * @tparam RealType the floating point number type
        * @param[in] l the parameter of the exponential distribution
        * @return a random number with an exponential distribution
        */
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType exp(RealType l)
        {
            auto exp_dist_l =
                std::exponential_distribution<RealType>(l);
            return exp_dist_l(m_rng);
        }

        /**
        * This function returns positive 64bit integers with
        * a Poisson's distribution (with parameter l)
        *
        * @tparam RealType the floating point number type of the parameter l
        * @param[in] l the parameter of the Poisson's distribution
        * @return a random integer with a Poisson's distribution
        */
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        uint64_t poisson(RealType l)
        {
            auto poisson_dist_l =
                std::poisson_distribution<uint64_t>(l);
            return poisson_dist_l(m_rng);
        }

        /**
        * This function returns floating point numbers with
        * a gaussian distribution (with given mean and deviation)
        *
        * @tparam RealType the floating point number type
        * @param[in] mean the mean of the gaussian distribution
        * @param[in] deviation the standard deviation of the gaussian distribution
        * @return a random number with a gaussian distribution
        */
        template<typename RealType>
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType gaussian(RealType mean, RealType deviation)
        {
            auto gaussian_dist_l =
                std::normal_distribution<RealType>(mean, deviation);
            return gaussian_dist_l(m_rng);
        }

    private:
        std::mt19937_64 m_rng;

    };

}
}
}

#endif //PICSAR_MULTIPHYSICS_RNG_WRAPPER
