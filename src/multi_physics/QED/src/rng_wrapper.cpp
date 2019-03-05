#include "rng_wrapper.h"

using namespace std;
using namespace picsar::multi_physics;

rng_wrapper::rng_wrapper(int64_t seed){
    rng.seed(seed);

    unf_dist_0_1 = bind(uniform_real_distribution<double>{0.0, 1.0}, ref(rng));
    exp_dist_l1 = bind(exponential_distribution<double>{}, ref(rng));
}

double rng_wrapper::get_unf_0_1(){
    return unf_dist_0_1();
}

double rng_wrapper::get_unf(double a, double b){
    auto unf_dist_a_b = bind(uniform_real_distribution<double>(a, b), ref(rng));
    return unf_dist_a_b();
}

double rng_wrapper::get_exp_l1(){
    return exp_dist_l1();
}
