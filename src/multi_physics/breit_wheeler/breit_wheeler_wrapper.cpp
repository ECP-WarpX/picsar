#include "breit_wheeler_wrapper.h"

using namespace multi_physics;

extern "C" void dummy_func_c_(double* a, double* b, double* res){
    *res = dummy_function(*a, *b);
}
