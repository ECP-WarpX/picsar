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
          //Minimum chi for photons to be considered by the engine
          _REAL chi_part_min =
            static_cast<_REAL>(__quantum_synchrotron_min_chi_part);

          _REAL chi_phot_min =
            static_cast<_REAL>(__quantum_synchrotron_min_chi_phot);

          _REAL chi_part_tdndt_min =
            static_cast<_REAL>(__quantum_synchrotron_min_tdndt_chi_part);
          _REAL chi_part_tdndt_max =
            static_cast<_REAL>(__quantum_synchrotron_max_tdndt_chi_part);
          size_t chi_part_tdndt_how_many =
            __quantum_synchrotron_how_many_tdndt_chi_part;
          tdndt_table_style tdndt_style =
            __quantum_synchrotron_dndt_table_style;
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
         quantum_synchrotron_engine(const quantum_synchrotron_engine& other);

         //Move constructor
         quantum_synchrotron_engine(quantum_synchrotron_engine&& other);

         //Getter & setter for lambda
         _REAL get_lambda() const;
         void set_lambda(_REAL lambda);

         //get a single optical depth
         PXRMP_FORCE_INLINE
         _REAL get_optical_depth();

         //Calculates the photon emission rate (Warning: no lookup tables)
         PXRMP_FORCE_INLINE
         _REAL compute_dN_dt(_REAL energy_part, _REAL chi_part) const;

         //Computes the lookup_table needed for dN/dt
         //(accepts pointer to ostream for diag)
         PXRMP_FORCE_INLINE
         void compute_dN_dt_lookup_table(std::ostream* stream = nullptr);

         //Interp the dN_dt from table
         PXRMP_FORCE_INLINE
         _REAL interp_dN_dt(_REAL energy_part, _REAL chi_part) const;


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
         const _REAL four = static_cast<_REAL>(4.0);
         const _REAL five = static_cast<_REAL>(5.0);

         //Internal function to compute KK function default
          //(accepts pointer to ostream for diag)
         void compute_KK_default_lookup_table(std::ostream* stream = nullptr);

         //Internal functions to perform calculations
         PXRMP_FORCE_INLINE
         _REAL compute_y(_REAL chi_phot, _REAL chi_ele) const;

         PXRMP_FORCE_INLINE
         _REAL compute_inner_integral(_REAL y) const;

         PXRMP_FORCE_INLINE
         _REAL compute_KK_integrand(_REAL chi_part, _REAL chi_phot) const;

         PXRMP_FORCE_INLINE
         _REAL compute_KK_function(_REAL chi_part) const;

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
quantum_synchrotron_engine(const quantum_synchrotron_engine& other):
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
    return rng.exp(one);
}

//Calculates the photon emission rate (Warning: no lookup tables)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_dN_dt(_REAL energy_part, _REAL chi_part) const
{
    if(energy_part == zero || chi_part == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__quantum_synchrotron_rate_coeff)*
        lambda*one/(energy_part*chi_part);

    return coeff*compute_KK_function(chi_part);
}

//Computes the lookup_table needed for dN/dt
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
void picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_dN_dt_lookup_table(std::ostream* stream)
{
    if(qs_ctrl.tdndt_style == tdnt_style_default)
        compute_KK_default_lookup_table(stream);

    //Other table styles are not currently implemented
}

//Interp the dN_dt from table
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
interp_dN_dt(_REAL energy_part, _REAL chi_part) const
{
    if(energy_part == zero || chi_part == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__quantum_synchrotron_rate_coeff)*
        lambda*one/(energy_part*chi_part);


    _REAL KK = zero;

    //If chi is out of table, use the first (or the last) value
    //in the table
    if(chi_part <= qs_ctrl.chi_part_tdndt_min){
        chi_part = qs_ctrl.chi_part_tdndt_min;
    }
    else if(chi_part >= qs_ctrl.chi_part_tdndt_max){
        chi_part = qs_ctrl.chi_part_tdndt_max ;
    }

    //Other table styles are not currently implemented
   if(qs_ctrl.tdndt_style == tdnt_style_default)
    KK =  exp(KKfunc_table.interp_linear(log(chi_part)));

    _REAL dndt = coeff * KK;

    return dndt;
}


//___________________PRIVATE FUNCTIONS_____________________________________

//Internal function to compute KK function default
 //(accepts pointer to ostream for diag)
template<typename _REAL, class _RNDWRAP>
void picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_KK_default_lookup_table(std::ostream* stream)
{
    //Prepare the KK_coords vector
    std::vector<_REAL> KK_coords = generate_log_spaced_vec(
     qs_ctrl.chi_part_tdndt_min, qs_ctrl.chi_part_tdndt_max,
    qs_ctrl.chi_part_tdndt_how_many);

    msg("Computing table for dNdt...\n", stream);
    //Do the hard work
    std::vector<_REAL> KK_vals{};
    std::transform(KK_coords.begin(), KK_coords.end(),
        std::back_inserter(KK_vals),
        [this](_REAL chi){return compute_KK_function(chi);});
    msg("...done!\n", stream);

    //The table will store the logarithms
    auto logfun = [](_REAL val){return log(val);};
    std::transform(KK_coords.begin(), KK_coords.end(), KK_coords.begin(), logfun);
    std::transform(KK_vals.begin(), KK_vals.end(), KK_vals.begin(), logfun);

    KKfunc_table = lookup_1d<_REAL>{KK_coords, KK_vals};
}

//Function to compute y (Warning: it doesn't chek if chi_ele != 0 or
//if chi_phot < chi_ele)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_y(_REAL chi_phot, _REAL chi_ele) const
{
    return chi_phot/(three*chi_ele*(chi_ele-chi_phot));
}

//Function to compute the inner integral of the QS photon emission rate
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_inner_integral(_REAL y) const
{
    auto func = [this](double y){
        return k_v(one/three, y);
    };
    return quad_a_inf<_REAL>(func, two*y);
}

template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_KK_integrand(_REAL chi_phot, _REAL chi_part) const
{
    if (chi_part == zero)
        return zero;

    _REAL y = compute_y(chi_phot, chi_part);

    _REAL inner = compute_inner_integral(y);

    _REAL part_2 = (two + three*chi_phot*y)*k_v(two/three, y);

    _REAL coeff = one/static_cast<_REAL>(pi*sqrt(three));

    return (inner - part_2)*coeff;
}

template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_KK_function(_REAL chi_part) const
{

    auto func = [chi_part, this](_REAL chi_phot){
        if(chi_part - chi_phot == zero || chi_phot == zero)
            return zero;
        else
            return compute_KK_integrand(chi_phot, chi_part);
    };

    return quad_a_b<_REAL>(func, qs_ctrl.chi_phot_min, chi_part);
}



#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
