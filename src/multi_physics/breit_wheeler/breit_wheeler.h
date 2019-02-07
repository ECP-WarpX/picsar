#ifndef __BREIT_WHEELER__
#define __BREIT_WHEELER__

#define __BUILD_WRAPPER__

namespace multi_physics{
    double dummy_function(double a, double b);

#ifdef __BUILD_WRAPPER__
    #include "breit_wheeler_wrapper.hpp"
#endif

}
#endif
