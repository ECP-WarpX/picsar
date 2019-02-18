#include "nonlin_breit_wheeler_engine.h"

using namespace std;
using namespace picsar::multi_physics;

nonlin_breit_wheeler_engine::nonlin_breit_wheeler_engine(int64_t seed, double _lambda):lambda(_lambda){
    rng.seed(seed);

    normalized_schwinger_field = calc_schwinger_given_lambda_SI(lambda);
    normalized_inv_schwinger_field = 1.0 / normalized_schwinger_field;
}
