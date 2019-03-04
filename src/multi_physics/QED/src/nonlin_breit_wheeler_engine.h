#ifndef __PMP_BREITWHEELER__
#define __PMP_BREITWHEELER__

//***********************WARNING! RNG from STL are NOT threadsafe! *************************

#include <cmath>
#include <cstdint>
#include <random>
#include <functional>
#include <fstream>

#include "commons.h"
#include "rng_wrapper.h"
#include "special_functions.h"
#include "quadrature.h"
#include "lookup_table.hpp"


namespace picsar{
  namespace multi_physics{

    typedef struct{
        double chi_phot_low;
        size_t chi_phot_how_many;
        double chi_phot_mul;
        double chi_ele_frac_min;
        size_t chi_ele_frac_how_many;
    } cumulative_distrib_params_list;

    typedef struct{
        double chi_phot_low;
        size_t chi_phot_how_many;
        double chi_phot_mul;
    } prod_rate_params_list;

    class nonlin_breit_wheeler_engine{
    public:
      nonlin_breit_wheeler_engine(int64_t seed, double lambda);

      //get a single optical depth
      double get_optical_depth();
      //initialize an optical depth vector
      void init_optical_depth_vector(std::vector<double>& opt_vec);

      double get_total_pair_production_rate(double gamma_phot, double chi_phot);

      void generate_tables(cumulative_distrib_params_list cum_params, prod_rate_params_list rate_params, std::ostream* diag = nullptr);
      void print_cumulative_distrib_pair_table(std::string file_name);

      bool has_lookup_tables();
    private:
        const double pair_prod_coeff =
            fine_structure_constant * electron_mass * speed_of_light * speed_of_light   / (reduced_plank * M_PI * sqrt(3.0));

        rng_wrapper rng;

        double lambda;
        double normalized_schwinger_field; //Normalized according to Smilei conventions
        double normalized_inv_schwinger_field;
        double rate_conversion_factor_SI_to_code;

        bool lookup_tables_flag = false;

        cumulative_distrib_params_list cumulative_distrib_params;
        lookup_table<2, double> cumulative_distrib_table;

        prod_rate_params_list prod_rate_params;
        lookup_table<1, double> T_table;

        static double compute_x(double chi_phot, double chi_ele);
        static double compute_inner_integral(double x);
        static double compute_TT_integrand(double chi_phot, double chi_ele);
        static double compute_TT_function(double chi_phot);
        static double compute_cumulative_distrib_pair(double chi_phot, double chi_ele);
        void generate_cumulative_distrib_pair_table(cumulative_distrib_params_list params, std::ostream* diag);
        void generate_TT_table(prod_rate_params_list params, std::ostream* diag);

        double compute_dN_dt(double gamma_phot, double chi_phot);

    };

  }
}

#endif
