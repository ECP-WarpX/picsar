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

      double unf(double a, double b); //Get rnd number uniformly distributed in [a,b)
      double expl1(); //Get rnd number with exponential distribution (lambda=1)

    private:
      std::mt19937_64 rng;

    };

  }
}


#endif
