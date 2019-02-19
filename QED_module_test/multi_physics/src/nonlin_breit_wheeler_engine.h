#ifndef __PMP_BREITWHEELER__
#define __PMP_BREITWHEELER__

//***********************WARNING! RNG from STL are NOT threadsafe! *************************

#include <cmath>
#include <cstdint>
#include <random>
#include <functional>

#include "commons.h"

namespace picsar{
  namespace multi_physics{

    class nonlin_breit_wheeler_engine{
    public:
      nonlin_breit_wheeler_engine(int64_t seed, double lambda);
      double calc_total_pair_production_rate(double chi_phot);

      bool has_lookup_tables();

      inline double get_unf_0_1(){
          return unf_dist_0_1();
      }

      inline double get_unf_eps_1(){
          return unf_dist_eps_1();
      }

      inline double get_optical_depth(){
          return -log(unf_dist_eps_1());
      }

    private:
        std::mt19937_64 rng;
        std::function<double()> unf_dist_0_1; //Uniform distribution between [0,1)
        std::function<double()> unf_dist_eps_1; //Uniform distrubtion between [e, 1), where e is numeric_limits<double>::min()

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
