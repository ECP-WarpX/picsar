#ifndef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
#define __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__

//This .hpp file contais the implementation of the
//nonlinear Breit Wheeler engine

#include <memory>

//Should be included by all the src files of the library
#include "qed_commons.h"

//Uses random numbers
#include "rng_wrapper.hpp"

//############################################### Declaration

namespace picsar{
  namespace multi_physics{

     //This class implements the Nonlinear Breit-Wheeler physical p
     //
     template<typename _REAL, class _RNDWRAP>
     class breit_wheeler_engine
     {
     public:
         //A unique pointer to a rng wrapper has to be passed to the
         //constructor. The constructor can accept a lambda parameter.
         //It is ignored if the SI units option is selected
         breit_wheeler_engine
         (_RNDWRAP&& rng,
         _REAL lambda = static_cast<_REAL>(1.0));

         //This explicit declaration makes sure that the move constructor
         //can be generated automatically
         breit_wheeler_engine(breit_wheeler_engine&& ) = default;

         //Getter & setter for lambda
         _REAL get_lambda() const;
         void set_lambda(_REAL lambda);

     private:
        _REAL lambda;
        _RNDWRAP rng;
     };

  }
}

//############################################### Implementation

template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
breit_wheeler_engine
(_RNDWRAP&& rng, _REAL lambda):
    lambda{lambda}, rng{std::move(rng)}
{
    //This enforces lambda=1 if SI units are used.
#ifdef PXRMP_WITH_SI_UNITS
    lambda = static_cast<_REAL>(1.0);
#endif
}


//Getter for lambda
template<typename _REAL, class _RNDWRAP>
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
get_lambda() const
{
    return lambda;
}

//Setter for lambda
template<typename _REAL, class _RNDWRAP>
void picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
set_lambda
#ifdef PXRMP_WITH_NORMALIZED_UNITS
(_REAL lambda)
{
    this->lambda = lambda;
}
#else
(_REAL){} //Do nothing
#endif

#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
