#ifndef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CTRL__
#define __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CTRL__

//Should be included by all the src files of the library
#include "qed_commons.h"

namespace picsar{
  namespace multi_physics{

      //This structure contains parameters which control how the BW engine
      //works
      template<typename _REAL>
      struct breit_wheeler_engine_ctrl{
           //Minimum chi_phot to consider
          _REAL chi_phot_min =
            static_cast<_REAL>(__breit_wheeler_min_chi_phot);

          _REAL chi_phot_tdndt_min =
            static_cast<_REAL>(__breit_wheeler_min_tdndt_chi_phot);
          _REAL chi_phot_tdndt_max =
            static_cast<_REAL>(__breit_wheeler_max_tdndt_chi_phot);
          size_t chi_phot_tdndt_how_many =
            __breit_wheeler_how_many_tdndt_chi_phot;

          _REAL chi_phot_tpair_min =
            static_cast<_REAL>(__breit_wheeler_min_tpair_chi_phot);
          _REAL chi_phot_tpair_max =
            static_cast<_REAL>(__breit_wheeler_max_tpair_chi_phot);
          size_t chi_phot_tpair_how_many =
            __breit_wheeler_how_many_tpair_chi_phot;
          size_t chi_frac_tpair_how_many =
            __breit_wheeler_chi_frac_tpair_how_many;
      };
  }
}

#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CTRL__
