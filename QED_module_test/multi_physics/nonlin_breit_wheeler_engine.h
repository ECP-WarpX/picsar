#ifndef __PMP_BREITWHEELER__
#define __PMP_BREITWHEELER__

#include <cmath>
#include <cstdint>
#include <random>

#include "commons.h"

namespace picsar{
  namespace multi_physics{

    class nonlin_breit_wheeler_engine{
    public:
      nonlin_breit_wheeler_engine(int64_t seed, double lambda);

    private:
        std::mt19937_64 rng;
        double lambda;
        double normalized_schwinger_field; //Normalized according to Smilei conventions
        double normalized_inv_schwinger_field;


    };

  }
}

#endif
