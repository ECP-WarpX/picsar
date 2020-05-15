#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar tables
#include "../../containers/picsar_tables.hpp"

#include "../../containers/picsar_array.hpp"

#include "../../containers/picsar_span.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    template<typename RealType>
    struct dndt_lookup_table_params{
        // dN/dt table:
        //__breit_wheeler_min_tdndt_chi_phot  is the inferior
        //limit of the total pair production rate lookup table.
        //If  __breit_wheeler_min_chi_phot < chi <__breit_wheeler_min_tdndt_chi_phot
        //BW process is taken into account, but the Erber approximation
        //rather than the lookup table is used.
        RealType chi_phot_min; //Min chi_phot
        RealType chi_phot_max; //Max chi_phot
        int chi_phot_how_many; //How many points
    };

    template<typename RealType, typename VectorType>
    class dndt_lookup_table{

        using view_type = dndt_lookup_table<RealType,
            containers::picsar_span<RealType>>;

        public:
            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(RealType where) const noexcept
            {
                return m_table.interp(where);
            }



            dndt_lookup_table_params<RealType, VectorType> m_ctrl;
            equispaced_1d_table<RealType, VectorType> m_table;
        };

}
}
}
}

#endif PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
