#include "nonlin_breit_wheeler_engine.h"

using namespace std;
using namespace picsar::multi_physics;

nonlin_breit_wheeler_engine::nonlin_breit_wheeler_engine(int64_t seed, double _lambda):
    rng{seed}, lambda{_lambda}{
    normalized_schwinger_field = calc_schwinger_given_lambda(lambda);
    normalized_inv_schwinger_field = 1.0 / normalized_schwinger_field;

    rate_conversion_factor_SI_to_code =  lambda/(2.0*speed_of_light*M_PI);
}

double nonlin_breit_wheeler_engine::get_optical_depth(){
    return rng.get_exp_l1();
}

void nonlin_breit_wheeler_engine::init_optical_depth_vector(std::vector<double>& opt_vec){
    for (auto& opt: opt_vec)
        opt = get_optical_depth();
}

bool nonlin_breit_wheeler_engine::has_lookup_tables(){
    return lookup_tables_flag;
}

double nonlin_breit_wheeler_engine::compute_inner_integral(double x){
    auto func = [](double s){return sqrt(s)*k_1_3(2.0*pow(s, 1.5)/3.0 ); };
    return quad_a_inf(func, x);
}

double nonlin_breit_wheeler_engine::compute_T_function(double chi_phot){
    if(chi_phot == 0)
        return 0.0;
    double coeff = 1./(M_PI * sqrt(3.0) * chi_phot * chi_phot);
    auto func = [chi_phot, this](double chi_ele){
        double X;
        if (chi_phot > chi_ele && chi_ele != 0)
            X = pow(chi_phot/(chi_ele*(chi_phot-chi_ele)), 2./3.);
        else
            X = BIG_POSITIVE_DOUBLE;
        double X_3_2 = pow(X, 3./2.);
        return compute_inner_integral(X)-(2.0-chi_phot*X_3_2)*k_2_3((2.0/3.0)*X_3_2);
    };
    return coeff*quad_a_b(func, 0.0, chi_phot);
}

// In CODE UNITS
double nonlin_breit_wheeler_engine::compute_dN_dt(double gamma_phot, double chi_phot){
     return rate_conversion_factor_SI_to_code*pair_prod_coeff*compute_T_function(chi_phot)*(chi_phot/gamma_phot);
 }

//Will use lookup tables soon
 double nonlin_breit_wheeler_engine::get_total_pair_production_rate(double gamma_phot, double chi_phot){
     return compute_dN_dt(gamma_phot, chi_phot);
 }
