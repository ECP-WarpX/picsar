#include "quadrature.h"

using namespace std;
using namespace picsar::multi_physics;

//If possible build with Boost library
#ifdef QUADRATURE_BUILD_WITH_BOOST

double picsar::multi_physics::quad_a_b(const std::function<double(double)>& func, double a, double b){
    return boost::math::quadrature::trapezoidal(func, a, b);
}

double picsar::multi_physics::quad_a_inf(const std::function<double(double)>& func, double a){ // I am really not happy with this...
    auto map_func = [a](double x) { return a - 1.0 +  1.0/ (1.0 - x); };
    auto inv_deriv_map_func = [](double x) { return 1/((1.0 - x)*(1.0 - x)); };
    auto mod_func = [&func, &map_func, &inv_deriv_map_func](double x) { return func(map_func(x))*inv_deriv_map_func(x); };
    return boost::math::quadrature::trapezoidal(mod_func, sqrt(numeric_limits<double>::epsilon()), 1.0 - sqrt(numeric_limits<double>::epsilon()));
}


//If not use fallback implementation (TO DO!)
#else

double picsar::multi_physics::quad_a_b(const std::function<double(double)>& func, double a, double b){
    //TO DO
    throw std::logic_error("Function not yet implemented");
    return 0;
}

double picsar::multi_physics::quad_a_inf(const std::function<double(double)>& func, double a){
    //TO DO
    throw std::logic_error("Function not yet implemented");
    return 0;
}

#endif
