#include "nonlin_breit_wheeler_engine.h"

using namespace std;
using namespace picsar::multi_physics;

nonlin_breit_wheeler_engine::nonlin_breit_wheeler_engine(int64_t seed, double _lambda):
    rng{seed}, lambda{_lambda}{
    normalized_schwinger_field = calc_schwinger_given_lambda(lambda);
    normalized_inv_schwinger_field = 1.0 / normalized_schwinger_field;

    rate_conversion_factor_SI_to_code =  lambda/(2.0*speed_of_light*M_PI);
}

double nonlin_breit_wheeler_engine::calc_total_pair_production_rate(double chi_phot){
    if(lookup_tables_flag)
        return calc_total_pair_production_rate_from_lookup(chi_phot);
    else
        return calc_total_pair_production_rate_from_scratch(chi_phot);
}

bool nonlin_breit_wheeler_engine::has_lookup_tables(){
    return lookup_tables_flag;
}

double nonlin_breit_wheeler_engine::compute_inner_integral(double x){
    auto func = [](double s){return sqrt(s)*k_1_3(2.0*pow(s, 1.5)/3.0 ); };
    return quad_a_inf(func, x);
}

double nonlin_breit_wheeler_engine::compute_d2N_dchi_dt(double gamma_phot, double chi_phot, double chi_ele){
    /*double chi_pos = chi_phot - chi_ele;
    double x = pow(chi_phot/(chi_pos*chi_ele), 2./3.);

    double inner = compute_inner_integral(x);
    double x_3_2 = pow(x, 3./2.);

    return rate_conversion_factor*pair_prod_coeff*(inner-(2.0-chi_phot*x_3_2)*k_2_3(2.0*x_3_2/3.0))/(gamma_phot*chi_phot);*/
    return 0.0;

}

double nonlin_breit_wheeler_engine::compute_T_function(double chi_phot){
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

// Now in SI units
 double nonlin_breit_wheeler_engine::compute_dN_dt(double gamma_phot, double chi_phot){
     return rate_conversion_factor_SI_to_code*pair_prod_coeff*compute_T_function(chi_phot)*(chi_phot/gamma_phot);
 }

double nonlin_breit_wheeler_engine::calc_total_pair_production_rate_from_scratch(double chi_phot){
    //TODO
    return 0.0;
}

double nonlin_breit_wheeler_engine::calc_total_pair_production_rate_from_lookup(double chi_phot){
    //TODO
    return 0.0;
}
