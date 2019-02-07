#ifndef __BREIT_WHEELER__
#define __BREIT_WHEELER__

#include "../common_types.h"

namespace multi_physics{
    double dummy_function(double a, double b);

    class breit_wheeler{
    public:
        static function_table_1d calculate_ritus_formula_table();
    };

}
#endif
