#include "nonlin_breit_wheeler_engine.h"

using namespace std;
using namespace picsar::multi_physics;

nonlin_breit_wheeler_engine::nonlin_breit_wheeler_engine(int64_t seed, double _lambda):lambda(_lambda){
    rng.seed(seed);

    unf_dist_0_1 = bind(uniform_real_distribution<double>{0.0, 1.0}, ref(rng));
    unf_dist_eps_1 = bind(uniform_real_distribution<double>{numeric_limits<double>::min(), 1.0}, ref(rng));

    normalized_schwinger_field = calc_schwinger_given_lambda_SI(lambda);
    normalized_inv_schwinger_field = 1.0 / normalized_schwinger_field;
}


double nonlin_breit_wheeler_engine::calc_total_pair_production_rate(double chi_phot){
    if(lookup_tables_flag)
        calc_total_pair_production_rate_from_lookup(chi_phot);
    else
        calc_total_pair_production_rate_from_scratch(chi_phot);
}

bool nonlin_breit_wheeler_engine::has_lookup_tables(){
    return lookup_tables_flag;
}


double nonlin_breit_wheeler_engine::calc_total_pair_production_rate_from_scratch(double chi_phot){
    //TODO
    return 0.0;
}

double nonlin_breit_wheeler_engine::calc_total_pair_production_rate_from_lookup(double chi_phot){
    //TODO
    return 0.0;
}
