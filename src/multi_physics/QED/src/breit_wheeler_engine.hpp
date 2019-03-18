#ifndef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
#define __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__

//This .hpp file contais the implementation of the
//nonlinear Breit Wheeler engine

#include <memory>

//Should be included by all the src files of the library
#include "qed_commons.h"

//Uses random numbers
#include "rng_wrapper.hpp"

//Uses special functions
#include "special_functions.hpp"

//Uses quadrature
#include "quadrature.hpp"

//############################################### Declaration

namespace picsar{
  namespace multi_physics{

      //Templates are used for the numerical type and for the
      //RNG wrapper
     template<typename _REAL, class _RNDWRAP>
     class breit_wheeler_engine
     {
     public:
         //A random number generatator has to be passed by move.
         //The RNG can be ANY object implementing the functions
         //_REAL unf (_REAL a, _REAL b)
         //and
         //_REAL exp (_REAL l)
         //The constructor can accept a lambda parameter.
         //It is ignored if the SI units option is selected
         breit_wheeler_engine
         (_RNDWRAP&& rng,
         _REAL lambda = static_cast<_REAL>(1.0));

         //This explicit declaration makes sure that the move constructor
         //can be generated automatically
          breit_wheeler_engine(breit_wheeler_engine& ):
              lambda(lambda), rng(rng)
              {};

         breit_wheeler_engine(breit_wheeler_engine&& ):
             lambda(std::move(lambda)), rng(std::move(rng))
             {};

         //Getter & setter for lambda
         _REAL get_lambda() const;
         void set_lambda(_REAL lambda);

         //get a single optical depth
         PXRMP_FORCE_INLINE
         _REAL get_optical_depth();

         //Calculates the pair production rate (Warning: from scratch!)
         PXRMP_FORCE_INLINE
         _REAL compute_dN_dt(_REAL energy_phot, _REAL chi_phot);

     private:
        _REAL lambda;

        //The only requrement for the RNG is to be able to provide unf(a,b) and
        //exp(l)
        _RNDWRAP rng;

        //Some handy constants
        const _REAL zero = static_cast<_REAL>(0.0);
        const _REAL one = static_cast<_REAL>(1.0);
        const _REAL two = static_cast<_REAL>(2.0);
        const _REAL three = static_cast<_REAL>(2.0);

        //Internal functions to perform calculations
        PXRMP_FORCE_INLINE
        _REAL compute_x(_REAL chi_phot, _REAL chi_ele) const;

        PXRMP_FORCE_INLINE
        _REAL compute_inner_integral(_REAL x) const;

        PXRMP_FORCE_INLINE
        _REAL compute_TT_integrand(_REAL chi_phot, _REAL chi_ele) const;

        PXRMP_FORCE_INLINE
        _REAL compute_TT_function(_REAL chi_phot) const;
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

//get a single optical depth (basically a single call to the internal RNG)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
get_optical_depth()
{
    return rng.exp(static_cast<_REAL>(1.0));
}

//Calculates the pair production rate (Warning: from scratch!)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_dN_dt(_REAL energy_phot, _REAL chi_phot)
{

    _REAL coeff = static_cast<_REAL>(__pair_prod_coeff)*
        lambda*(one/( chi_phot * energy_phot));
    return coeff*compute_TT_function(chi_phot);
}

//Function to compute X (Warning: it doen't chek if chi_ele != 0 or
//if chi_phot > chi_ele)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_x(_REAL chi_phot, _REAL chi_ele) const
{
    return pow(chi_phot/(chi_ele*(chi_phot-chi_ele)), two/three);
}

//Function to compute the inner integral of the pair production rate
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_inner_integral(_REAL x) const
{
    auto func = [this](double s){
        return sqrt(s)*k_v(two/three, two*pow(s, three/two)/three);
    };
    return quad_a_inf<_REAL>(func, x);
}

//Calculations of other parts of the pair production rate [I]
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_TT_integrand(_REAL chi_phot, _REAL chi_ele) const
{
    _REAL xx = compute_x(chi_phot, chi_ele);
    _REAL xx_3_2 = pow(xx, three/two);
    return compute_inner_integral(xx)-(two-chi_phot*xx_3_2)
        *k_v(two/three, (two/three)*xx_3_2);
}

//Calculations of other parts of the pair production rate [II]
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_TT_function(_REAL chi_phot) const
{
    if(chi_phot == zero)
        return zero;

    auto func = [chi_phot, this](_REAL chi_ele){
        if(chi_ele - chi_phot == zero)
            return zero;
        else
            return compute_TT_integrand(chi_phot, chi_ele);
    };

    return quad_a_b<_REAL>(func, zero, chi_phot);
}

#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
