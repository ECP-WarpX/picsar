#ifndef __PMP_BREITWHEELER__
#define __PMP_BREITWHEELER__

//***********************WARNING! RNG from STL are NOT threadsafe! *************************

#include <cmath>
#include <cstdint>
#include <random>
#include <functional>

#include "commons.h"
#include "rng_wrapper.h"

namespace picsar{
  namespace multi_physics{

    class nonlin_breit_wheeler_engine{
    public:
      nonlin_breit_wheeler_engine(int64_t seed, double lambda);
      double calc_total_pair_production_rate(double chi_phot);

      bool has_lookup_tables();

    private:
        rng_wrapper rng;

        double lambda;
        double normalized_schwinger_field; //Normalized according to Smilei conventions
        double normalized_inv_schwinger_field;

        bool lookup_tables_flag = false;

        double calc_total_pair_production_rate_from_scratch(double chi_phot);

        double calc_total_pair_production_rate_from_lookup(double chi_phot);


    };

  }
}

#endif
