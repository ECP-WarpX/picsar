#ifndef __PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CTRL__
#define __PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CTRL__

//Should be included by all the src files of the library
#include "qed_commons.h"


namespace picsar{
  namespace multi_physics{

      //This structure contains parameters which control how the QS engine
      //works
      template<typename _REAL>
      struct quantum_synchrotron_engine_ctrl{
          //Minimum chi for particles to be considered by the engine
          _REAL chi_part_min =
            static_cast<_REAL>(__quantum_synchrotron_min_chi_part);

          _REAL chi_part_tdndt_min =
            static_cast<_REAL>(__quantum_synchrotron_min_tdndt_chi_part);
          _REAL chi_part_tdndt_max =
            static_cast<_REAL>(__quantum_synchrotron_max_tdndt_chi_part);
          size_t chi_part_tdndt_how_many =
            __quantum_synchrotron_how_many_tdndt_chi_part;

          _REAL chi_part_tem_min =
            static_cast<_REAL>(__quantum_synchrotron_min_tem_chi_part);
          _REAL chi_part_tem_max =
            static_cast<_REAL>(__quantum_synchrotron_max_tem_chi_part);
          size_t chi_part_tem_how_many =
            __quantum_synchrotron_how_many_tem_chi_part;
          size_t prob_tem_how_many =
            __quantum_synchrotron_prob_tem_how_many;
      };
  }
}
#endif //__PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CTRL__