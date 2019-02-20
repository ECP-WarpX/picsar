#ifndef __PMP_RNG_WRAPPER__
#define __PMP_RNG_WRAPPER__

#include <random>
#include <functional>
#include <cstdint>

namespace picsar{
  namespace multi_physics{

    class rng_wrapper{
    public:
      rng_wrapper(int64_t seed);


      double get_unf_0_1(); //Get rnd number uniformly distributed in [0,1)
      double get_exp_l1(); //Get rnd number with exponential distribution (lambda=1)

    private:
      std::mt19937_64 rng;

      std::function<double()> unf_dist_0_1; //Uniform distribution between [0,1)
      std::function<double()> exp_dist_l1; //Exponential distribution, with lambda=1

    };

  }
}


#endif
