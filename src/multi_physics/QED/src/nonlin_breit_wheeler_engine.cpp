#include "nonlin_breit_wheeler_engine.h"

using namespace std;
using namespace std::placeholders;
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

double nonlin_breit_wheeler_engine::compute_x(double chi_phot, double chi_ele){
    if (chi_phot > chi_ele && chi_ele != 0)
        return pow(chi_phot/(chi_ele*(chi_phot-chi_ele)), 2./3.);
    else
        return BIG_POSITIVE_DOUBLE;
}

double nonlin_breit_wheeler_engine::compute_inner_integral(double x){
    auto func = [](double s){
        return sqrt(s)*k_1_3(2.0*pow(s, 1.5)/3.0 );
    };
    return quad_a_inf(func, x);
}

double nonlin_breit_wheeler_engine::compute_TT_integrand(double chi_phot, double chi_ele){
    double X = compute_x(chi_phot, chi_ele);
    double X_3_2 = pow(X, 3./2.);
    return compute_inner_integral(X)-(2.0-chi_phot*X_3_2)*k_2_3((2.0/3.0)*X_3_2);
}

double nonlin_breit_wheeler_engine::compute_TT_function(double chi_phot){
    if(chi_phot == 0)
        return 0.0;
    return quad_a_b(std::bind(compute_TT_integrand, chi_phot, _1), 0.0, chi_phot);
}

void nonlin_breit_wheeler_engine::generate_tables(cumulative_distrib_params_list cum_params, prod_rate_params_list prod_params, std::ostream* diag){
    generate_cumulative_distrib_pair_table(cum_params, diag);

    generate_pair_prod_table(prod_params, diag);

    lookup_tables_flag = true;
}

void nonlin_breit_wheeler_engine::generate_cumulative_distrib_pair_table(cumulative_distrib_params_list params, std::ostream* diag){
    message("    Generation of cumulative_distrib_table: START", diag);
    double chi_phot = params.chi_phot_low;

    std::vector<double> chi_phot_v(params.chi_phot_how_many);
    std::vector<double> chi_ele_frac_v(params.chi_ele_frac_how_many);

    generate(chi_phot_v.begin(), chi_phot_v.end(),
        [&chi_phot, params](){ double elem = chi_phot; chi_phot*=params.chi_phot_mul; return elem;}
    );

    size_t imax = chi_ele_frac_v.size()/2;
    double alpha = -atanh(params.chi_ele_frac_min*2.0 - 1.0)/imax;
    for(auto it = chi_ele_frac_v.begin(); it !=  chi_ele_frac_v.end(); it++){
        int iii = std::distance(chi_ele_frac_v.begin(), it) - imax;
        double ddd = iii*alpha;
        *it = (tanh(ddd) + 1.0)*0.5;
    }

    cumulative_distrib_params = params;

    cumulative_distrib_table = lookup_table<2, double>{chi_phot_v, chi_ele_frac_v};
    for(auto cpp: chi_phot_v){
        for(auto eef: chi_ele_frac_v){
            cumulative_distrib_table.fill_at(std::array<double,2>{cpp, eef}, compute_cumulative_distrib_pair(cpp, cpp*eef));
        }
        message("    Generation of cumulative_distrib_table: chi_phot = " + std::to_string(cpp), diag);
    }
    message("    Generation of cumulative_distrib_table: END", diag);
}
void nonlin_breit_wheeler_engine::generate_pair_prod_table(prod_rate_params_list params, std::ostream* diag){
    message("    Generation of pair_prod_table: START", diag);
    // TO DO
    message("    Generation of pair_prod_table: END", diag);
}

void nonlin_breit_wheeler_engine::print_cumulative_distrib_pair_table(std::string file_name){
    cumulative_distrib_table.print_on_disk(file_name);
}

// In CODE UNITS
double nonlin_breit_wheeler_engine::compute_dN_dt(double gamma_phot, double chi_phot){
    return rate_conversion_factor_SI_to_code*pair_prod_coeff*(1.0/( chi_phot * chi_phot))*compute_TT_function(chi_phot)*(chi_phot/gamma_phot);
 }

 double nonlin_breit_wheeler_engine::compute_cumulative_distrib_pair(double chi_phot, double chi_ele){
    double num = quad_a_b(std::bind(compute_TT_integrand, chi_phot, _1), 0.0, chi_ele);
    return num/compute_TT_function(chi_phot) ;
 }

//Will use lookup tables soon
 double nonlin_breit_wheeler_engine::get_total_pair_production_rate(double gamma_phot, double chi_phot){
     return compute_dN_dt(gamma_phot, chi_phot);
 }
