#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES_GENERATOR
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES_GENERATOR

//Should be included by all the src files of the library
#include "../../qed_commons.h"

#include "breit_wheeler_engine_tables.hpp"
#include "breit_wheeler_engine_tabulated_functions.hpp"

#include <omp.h>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    template<
        typename RealType,
        typename VectorType>
    void
    dndt_lookup_table<RealType, VectorType>::generate()
    {
        
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES_GENERATOR