#include "rng_wrapper.h"

using namespace std;
using namespace picsar::multi_physics;

rng_wrapper::rng_wrapper(int64_t seed){
    rng.seed(seed);
}

double rng_wrapper::unf(double a, double b){
    auto unf_dist_a_b = bind(uniform_real_distribution<double>(a, b), ref(rng));
    return unf_dist_a_b();
}

double rng_wrapper::expl1(){
    auto exp_dist_l1 = bind(exponential_distribution<double>{}, ref(rng));
    return exp_dist_l1();
}
