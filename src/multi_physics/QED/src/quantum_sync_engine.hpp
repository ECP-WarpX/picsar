#ifndef __PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE__
#define __PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE__

//This .hpp file contais the implementation of the
//quantum synchrotron emission engine

//Should be included by all the src files of the library
#include "qed_commons.h"

//Uses random numbers
#include "rng_wrapper.hpp"

//Uses special functions
#include "special_functions.hpp"

//Uses vec functions
#include "vec_functions.hpp"

//Uses quadrature
#include "quadrature.hpp"

//Uses chi functions
#include "chi_functions.hpp"

//Uses lookup tables
#include "lookup_tables.hpp"

//Uses messages
#include "msg.hpp"

//############################################### Declaration

namespace picsar{
  namespace multi_physics{

      //This structure contains parameters which control how the QS engine
      //works
      template<typename _REAL>
      struct quantum_synchrotron_engine_ctrl{
           //Minimum chi_phot to consider
          _REAL nu_min =static_cast<_REAL>(__qs_min_nu);
      };

      //Templates are used for the numerical type and for the
      //RNG wrapper
     template<typename _REAL, class _RNDWRAP>
     class quantum_synchrotron_engine
     {
     public:
         //A random number generatator has to be passed by move.
         //The RNG can be ANY object implementing the functions
         //_REAL unf (_REAL a, _REAL b)
         //and
         //_REAL exp (_REAL l)
         //The constructor can accept a lambda parameter.
         //It is ignored if the SI units option is selected
         //The constructor can accept also a quantum_synchrotron_engine_ctrl
         //struct, which controls how the engine work
         quantum_synchrotron_engine
         (_RNDWRAP&& rng,
         _REAL lambda = static_cast<_REAL>(1.0),
         quantum_synchrotron_engine_ctrl<_REAL> bw_ctrl =
         quantum_synchrotron_engine_ctrl<_REAL>());

         //Copy constructor
         quantum_synchrotron_engine(quantum_synchrotron_engine& other);

         //Move constructor
         quantum_synchrotron_engine(quantum_synchrotron_engine&& other);

         //Getter & setter for lambda
         _REAL get_lambda() const;
         void set_lambda(_REAL lambda);

         //get a single optical depth
         PXRMP_FORCE_INLINE
         _REAL get_optical_depth();
     private:
         _REAL lambda;

         //The only requrement for the RNG is to be able to provide unf(a,b) and
         //exp(l)
         _RNDWRAP rng;

         //Parameters which control how the engine works
         const quantum_synchrotron_engine_ctrl<_REAL> qs_ctrl;

         //lookup table for the KK function
         lookup_1d<_REAL> KKfunc_table;

         //lookup table for the cumulativie distribution table
         lookup_2d<_REAL> cum_distrib_table;
         //Auxiliary table for coordinate interpolation
         lookup_1d<_REAL> aux_table;

         //Some handy constants
         const _REAL zero = static_cast<_REAL>(0.0);
         const _REAL one = static_cast<_REAL>(1.0);
         const _REAL two = static_cast<_REAL>(2.0);
         const _REAL three = static_cast<_REAL>(3.0);
     };
  }
}

//############################################### Implementation


template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
quantum_synchrotron_engine
(_RNDWRAP&& rng, _REAL lambda, quantum_synchrotron_engine_ctrl<_REAL> qs_ctrl):
    lambda{lambda}, rng{std::move(rng)}, qs_ctrl{qs_ctrl}
{
    //This enforces lambda=1 if SI units are used.
#ifdef PXRMP_WITH_SI_UNITS
    lambda = static_cast<_REAL>(1.0);
#endif
}

//Copy constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
quantum_synchrotron_engine(quantum_synchrotron_engine& other):
    lambda(other.lambda), rng(other.rng), qs_ctrl(other.qs_ctrl),
    KKfunc_table(other.KKfunc_table), cum_distrib_table(other.cum_distrib_table),
    aux_table(other.aux_table)
    {}

//Move constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
quantum_synchrotron_engine(quantum_synchrotron_engine&& other):
    lambda(std::move(other.lambda)), rng(std::move(other.rng)),
    qs_ctrl(std::move(other.qs_ctrl)),
    KKfunc_table(std::move(other.KKfunc_table)),
    cum_distrib_table(std::move(other.cum_distrib_table)),
    aux_table(std::move(other.aux_table))
    {}


//Getter for lambda
template<typename _REAL, class _RNDWRAP>
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
get_lambda() const
{
    return lambda;
}

//Setter for lambda
template<typename _REAL, class _RNDWRAP>
void picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
set_lambda
#ifdef PXRMP_WITH_NORMALIZED_UNITS
(_REAL lambda)
{
    this->lambda = lambda;
}
#else
(_REAL){} //Do nothing
#endif

//get a single optical depth (basically a single call to the internal RNG)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
get_optical_depth()
{
    return rng.exp(static_cast<_REAL>(1.0));
}

//___________________PRIVATE FUNCTIONS_____________________________________

#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
