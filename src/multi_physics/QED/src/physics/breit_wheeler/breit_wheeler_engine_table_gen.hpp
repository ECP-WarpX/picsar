#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLE_GEN
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLE_GEN

#include <omp.h>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses physical constants
#include "breit_wheeler_engine_tables.hpp"

//Uses physical constants
#include "breit_wheeler_engine_tabulated_functions.hpp"

#include "breit_wheeler_engine_core.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    template<typename RealType, typename VectorType>
    dndt_lookup_table<RealType, VectorType>
    fill_dndt_lookup_table(dndt_lookup_table_params<RealType> params)
    {
    }
}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLE_GEN
