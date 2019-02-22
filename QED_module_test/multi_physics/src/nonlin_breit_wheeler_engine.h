#ifndef __PMP_BREITWHEELER__
#define __PMP_BREITWHEELER__

//***********************WARNING! RNG from STL are NOT threadsafe! *************************

#include <cmath>
#include <cstdint>
#include <random>
#include <functional>

#include "commons.h"
#include "rng_wrapper.h"
#include "special_functions.h"
#include "quadrature.h"


namespace picsar{
  namespace multi_physics{

    class nonlin_breit_wheeler_engine{
    public:
      nonlin_breit_wheeler_engine(int64_t seed, double lambda);
      double calc_total_pair_production_rate(double chi_phot);

      bool has_lookup_tables();

      double compute_d2N_dchi_dt(double gamma_phot, double chi_phot, double chi_ele);
      double compute_dN_dt(double gamma_phot, double chi_phot);

    private:
        const double pair_prod_coeff =
            fine_structure_constant * electron_mass * speed_of_light * speed_of_light   / reduced_plank;

        rng_wrapper rng;

        double lambda;
        double normalized_schwinger_field; //Normalized according to Smilei conventions
        double normalized_inv_schwinger_field;
        double rate_conversion_factor_SI_to_code;


        bool lookup_tables_flag = false;

        double compute_inner_integral(double x);

        double compute_T_function(double chi_phot);

        double calc_total_pair_production_rate_from_scratch(double chi_phot);

        double calc_total_pair_production_rate_from_lookup(double chi_phot);


    };

  }
}

#endif
