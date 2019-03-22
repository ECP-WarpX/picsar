#ifndef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
#define __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__

//This .hpp file contais the implementation of the
//nonlinear Breit Wheeler engine

#include <memory>
#include <utility>
#include<iostream>

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

//Uses utilities
#include "utilities.hpp"

//############################################### Declaration

namespace picsar{
  namespace multi_physics{

      //This enumerator contains all the possible table styles

      //This structure contains parameters which control how the BW engine
      //works
      template<typename _REAL>
      struct breit_wheeler_engine_ctrl{
           //Minimum chi_phot to consider
          _REAL chi_phot_min =static_cast<_REAL>(__bw_min_chi_phot);

          _REAL chi_phot_tdndt_min =static_cast<_REAL>(__bw_min_tdndt_chi_phot);
          _REAL chi_phot_tdndt_max =static_cast<_REAL>(__bw_max_tdndt_chi_phot);
          size_t chi_phot_tdndt_how_many =__bw_how_many_tdndt_chi_phot;
          tdndt_table_style tdndt_style = __bw_dndt_table_style;

          _REAL chi_phot_tpair_min =static_cast<_REAL>(__bw_min_tpair_chi_phot);
          _REAL chi_phot_tpair_max =static_cast<_REAL>(__bw_max_tpair_chi_phot);
          size_t chi_phot_tpair_how_many =__bw_how_many_tpair_chi_phot;

          size_t chi_frac_tpair_how_many = __chi_frac_tpair_how_many;

          tpair_table_style tpair_style = __bw_pair_table_style;
      };

      //This struct is used for post-run check of the suitability of the
      //pair table.
      struct breit_wheeler_error_info{
          size_t pair_prod_chi_too_big_how_many = 0;
          size_t pair_prod_chi_too_small_how_many = 0;
      };

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
         //The constructor can accept also a breit_wheeler_engine_ctrl
         //struct, which controls how the engine work
         breit_wheeler_engine
         (_RNDWRAP&& rng,
         _REAL lambda = static_cast<_REAL>(1.0),
         breit_wheeler_engine_ctrl<_REAL> bw_ctrl =
         breit_wheeler_engine_ctrl<_REAL>());

         //Copy constructor
         breit_wheeler_engine(breit_wheeler_engine& other);

         //Move constructor
         breit_wheeler_engine(breit_wheeler_engine&& other);

         //Getter & setter for lambda
         _REAL get_lambda() const;
         void set_lambda(_REAL lambda);

         //get a single optical depth
         PXRMP_FORCE_INLINE
         _REAL get_optical_depth();

         //Calculates the pair production rate (Warning: no lookup tables)
         PXRMP_FORCE_INLINE
         _REAL compute_dN_dt(_REAL energy_phot, _REAL chi_phot) const;

         //Computes the lookup_table needed for dN/dt
         //(accepts pointer to ostream for diag)
         PXRMP_FORCE_INLINE
         void compute_dN_dt_lookup_table(std::ostream* stream = nullptr);

         //Interp the dN_dt from table
         PXRMP_FORCE_INLINE
         _REAL interp_dN_dt(_REAL energy_phot, _REAL chi_phot) const;

         //This function evolves the optical depth for a particle and
         //checks if it goes to zero. If it doesn't the output is false,0.
         //On the contrary, if it goes to zero the output is true, dt_prod,
         //where dt_prod (<= dt) is the time at which the event occurs.
         PXRMP_FORCE_INLINE
         std::pair<bool, _REAL> evolve_opt_depth_and_determine_event(
         _REAL px, _REAL py, _REAL pz,
         _REAL ex, _REAL ey, _REAL ez,
         _REAL bx, _REAL by, _REAL bz,
         _REAL dt, _REAL& opt_depth) const;

         //Computes the cumulative pair production rate given
         //chi_phot and chi_ele
         PXRMP_FORCE_INLINE
         _REAL compute_cumulative_pair (_REAL chi_phot, _REAL chi_ele) const;

         //Computes the cumulative pair production rate lookup table
         void compute_cumulative_pair_table (std::ostream* stream = nullptr);

         //TO WRITE....
         PXRMP_FORCE_INLINE
         std::array<std::vector<std::pair<vec3<_REAL>, _REAL>>,2>
          generate_breit_wheeler_pairs(
         _REAL px, _REAL py, _REAL pz,
         _REAL ex, _REAL ey, _REAL ez,
         _REAL bx, _REAL by, _REAL bz,
         _REAL weight, size_t sampling) ;

         //get a copy of the struct storing info on the out of table errors
         breit_wheeler_error_info get_error_info();

     private:
        _REAL lambda;

        //The only requrement for the RNG is to be able to provide unf(a,b) and
        //exp(l)
        _RNDWRAP rng;

        //Parameters which control how the engine works
        const breit_wheeler_engine_ctrl<_REAL> bw_ctrl;

        //lookup table for the TT function
        lookup_1d<_REAL> TTfunc_table;

        //lookup table for the cumulativie distribution table
        lookup_2d<_REAL> cum_distrib_table;
        //Auxiliary table for coordinate interpolation
        lookup_1d<_REAL> aux_table;

        //This struct stores information relevant for post-simultation checks
        breit_wheeler_error_info err_info;

        //Some handy constants
        const _REAL zero = static_cast<_REAL>(0.0);
        const _REAL one = static_cast<_REAL>(1.0);
        const _REAL two = static_cast<_REAL>(2.0);
        const _REAL three = static_cast<_REAL>(3.0);

        //Internal function to compute TT function default
         //(accepts pointer to ostream for diag)
        void compute_TT_default_lookup_table(std::ostream* stream = nullptr);

        //Internal function to compute pair production table
        //(accepts pointer to ostream for diag)
        void compute_pair_default_lookup_table(std::ostream* stream = nullptr);


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
(_RNDWRAP&& rng, _REAL lambda, breit_wheeler_engine_ctrl<_REAL> bw_ctrl):
    lambda{lambda}, rng{std::move(rng)}, bw_ctrl{bw_ctrl}
{
    //This enforces lambda=1 if SI units are used.
#ifdef PXRMP_WITH_SI_UNITS
    lambda = static_cast<_REAL>(1.0);
#endif
}

//Copy constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
breit_wheeler_engine(breit_wheeler_engine& other):
    lambda(other.lambda), rng(other.rng), bw_ctrl(other.bw_ctrl),
    TTfunc_table(other.TTfunc_table), cum_distrib_table(other.cum_distrib_table),
    aux_table(other.aux_table),
    err_info(other.err_info)
    {}

//Move constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
breit_wheeler_engine(breit_wheeler_engine&& other):
    lambda(std::move(other.lambda)), rng(std::move(other.rng)),
    bw_ctrl(std::move(other.bw_ctrl)),
    TTfunc_table(std::move(other.TTfunc_table)),
    cum_distrib_table(std::move(other.cum_distrib_table)),
    aux_table(std::move(other.aux_table)),
    err_info(other.err_info)
    {}


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

//Calculates the pair production rate (Warning: no tables are used)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_dN_dt(_REAL energy_phot, _REAL chi_phot) const
{
    if(energy_phot == zero || chi_phot == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__pair_prod_coeff)*
        lambda*(one/( chi_phot * energy_phot));
    return coeff*compute_TT_function(chi_phot);
}

//Computes the lookup_table needed for dN/dt
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
void picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_dN_dt_lookup_table(std::ostream* stream)
{
    if(bw_ctrl.tdndt_style == tdnt_style_default)
        compute_TT_default_lookup_table(stream);

    //Other table styles are not currently implemented
}

//Interp the dN_dt from table
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
interp_dN_dt(_REAL energy_phot, _REAL chi_phot) const
{
    if(energy_phot == zero || chi_phot == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__pair_prod_coeff)*
        lambda*(one/( chi_phot * energy_phot));


        _REAL TT = zero;
    //Use approximate analytical expression if chi < chi_phot_tdndt_min
    //or chi > chi_phot_tdndt_min
    if(chi_phot <= bw_ctrl.chi_phot_tdndt_min ||
        chi_phot >= bw_ctrl.chi_phot_tdndt_max){
        _REAL a = static_cast<_REAL>(__erber_Tfunc_asynt_a);
        _REAL b = static_cast<_REAL>(__erber_Tfunc_asynt_b);
        TT = a*chi_phot*pow(k_v(one/three, b/chi_phot),two);
    }
    //otherwise use lookup tables
    else{
        //Other table styles are not currently implemented
        if(bw_ctrl.tdndt_style == tdnt_style_default)
            TT =  exp(TTfunc_table.interp(log(chi_phot)));
    }
    //**end

    _REAL dndt = coeff * TT;

    return dndt;
}


//This function evolves the optical depth for a particle and
//checks if it goes to zero. If it doesn't the output is false,0.
//On the contrary, if it goes to zero the output is true, dt_em,
//where dt_em (<= dt) is the time at which the event occurs.
//Warning! For now it does not use tables!
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
std::pair<bool, _REAL>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
evolve_opt_depth_and_determine_event(
_REAL px, _REAL py, _REAL pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dt, _REAL& opt_depth) const
{
    _REAL energy = norm<_REAL>(vec3<_REAL> {px, py, pz})*__c;
    _REAL chi = chi_photon(px, py, pz, ex, ey, ez, bx, by, bz, lambda);

    auto false_zero_pair = std::make_pair(false, zero);

    //Do NOT evolve opt_depth if the chi parameter is less then threshold
    //or if the photon energy is not high enough to generate a pair
    if(chi <= bw_ctrl.chi_phot_min ||
       energy < two*static_cast<_REAL>(__emass*__c*__c))
            return false_zero_pair;

    //**Compute dndt
    _REAL dndt;
    //Uses table if available
    if(TTfunc_table.is_init()){
        dndt = interp_dN_dt(energy, chi);
    }
    //If not it computes dndt
    else{
        err("dndt lookup table not initialized!\n");
        dndt = compute_dN_dt(energy, chi);
    }

    opt_depth -= dndt*dt;

    if(opt_depth > zero){
        return false_zero_pair;
    }
    else{
        //Calculates the time at which pair prod takes place
        _REAL dt_prod = opt_depth/dndt + dt;
        return std::make_pair(true, dt_prod);
    }
}


template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_cumulative_pair_table (std::ostream* stream)
{
    if(bw_ctrl.tpair_style == tpair_style_default)
        compute_pair_default_lookup_table(stream);

    //Other table styles are not currently implemented
}

//TO WRITE ...
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
std::array<std::vector<std::pair<picsar::multi_physics::vec3<_REAL>, _REAL>>,2>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
generate_breit_wheeler_pairs(
_REAL px, _REAL py, _REAL pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL weight, size_t sampling)
{
    std::vector<std::pair<vec3<_REAL>, _REAL>> electrons(sampling);
    std::vector<std::pair<vec3<_REAL>, _REAL>> positrons(sampling);

    _REAL chi_phot = chi_photon(px, py, pz, ex, ey, ez, bx, by, bz, lambda);

    //TO DO : casi speciali al limite della tavola

    size_t i = 0;
    for (auto frac: aux_table.ref_coords()){
        aux_table.ref_data()[i] = cum_distrib_table.interp(chi_phot, frac); //Definitely not the smartest thing to do....
        i++;
    }

    _REAL new_weight = weight/sampling;

    _REAL me_c = static_cast<_REAL>(__emass*__c);

    vec3<_REAL> p_phot{px, py, pz};
    _REAL norm_phot = norm(p_phot);
    vec3<_REAL> n_phot = p_phot/norm_phot;
    _REAL gamma_phot = norm_phot/me_c;

    for(size_t s = 0; s < sampling; s++){
        _REAL prob = rng.unf(zero, one);
        bool invert = false;
        if(prob > one/two){
            prob = one - prob;
            invert = true;
        }

        _REAL chi_ele = aux_table.interp(prob);
        _REAL chi_pos = chi_phot - chi_ele;

        if(invert)
            std::swap(chi_ele, chi_pos);

        _REAL coeff =  (gamma_phot - two)/chi_phot;

        _REAL cc_ele = sqrt((one+chi_ele*coeff)*(one+chi_ele*coeff)-one)*me_c;
        _REAL cc_pos = sqrt((one+chi_pos*coeff)*(one+chi_pos*coeff)-one)*me_c;


        vec3<_REAL> p_ele = cc_ele*n_phot;
        vec3<_REAL> p_pos = cc_pos*n_phot;

        electrons[s] = std::make_pair(p_ele, new_weight);
        positrons[s] = std::make_pair(p_pos, new_weight);
    }
    return {electrons, positrons};
}

//get a copy of the struct storing info on the out of table errors
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_error_info
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::get_error_info()
{
    return err_info;
}

//___________________PRIVATE FUNCTIONS_____________________________________

//Internal function to compute T function _default_
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_TT_default_lookup_table(std::ostream* stream)
{
    //Prepare the TT_coords vector
    std::vector<_REAL> TT_coords = generate_log_spaced_vec(
     bw_ctrl.chi_phot_tdndt_min, bw_ctrl.chi_phot_tdndt_max,
    bw_ctrl.chi_phot_tdndt_how_many);

    msg("Computing table for dNdt...\n", stream);
    //Do the hard work
    std::vector<_REAL> TT_vals{};
    std::transform(TT_coords.begin(), TT_coords.end(),
        std::back_inserter(TT_vals),
        [this](_REAL chi){return compute_TT_function(chi);});
    msg("...done!\n", stream);

    //The table will store the logarithms
    auto logfun = [](_REAL val){return log(val);};
    std::transform(TT_coords.begin(), TT_coords.end(), TT_coords.begin(), logfun);
    std::transform(TT_vals.begin(), TT_vals.end(), TT_vals.begin(), logfun);

    TTfunc_table = lookup_1d<_REAL>{TT_coords, TT_vals,
        lookup_1d<_REAL>::linear_interpolation};
}

//Internal function to compute pair production _default_ table
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_pair_default_lookup_table(std::ostream* stream)
{
    //Prepare the chi_coords vector
    std::vector<_REAL> chi_coords = generate_log_spaced_vec(
     bw_ctrl.chi_phot_tpair_min, bw_ctrl.chi_phot_tpair_max,
    bw_ctrl.chi_phot_tpair_how_many);

    std::vector<_REAL> frac_coords = generate_lin_spaced_vec(zero,one/two,
    bw_ctrl.chi_frac_tpair_how_many);

    std::vector<_REAL> pair_vals{};

    msg("Computing table for pair production...\n", stream);

    for(auto chi_phot: chi_coords){
        pair_vals.push_back(zero);
        msg("chi_phot: " + std::to_string(chi_phot) + " \n", stream);
        for(size_t i = 1; i < frac_coords.size() - 1; i++){
            _REAL temp = compute_cumulative_pair(
                chi_phot, chi_phot*frac_coords[i]);
            pair_vals.push_back(temp);
        }
        pair_vals.push_back(one/two); //The function is symmetric
    }
    msg("...done!\n", stream);

    cum_distrib_table = lookup_2d<_REAL>{
        std::array<std::vector<_REAL>,2>{chi_coords, frac_coords}, pair_vals,
        lookup_2d<_REAL>::linear_interpolation, lookup_2d<_REAL>::row_major};


    //Initialize the auxiliary table
    aux_table = lookup_1d<_REAL>{cum_distrib_table.get_coords()[1],
    std::vector<_REAL>(bw_ctrl.chi_frac_tpair_how_many),
    lookup_1d<_REAL>::linear_interpolation};
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
        return sqrt(s)*k_v(one/three, two*pow(s, three/two)/three);
    };
    return quad_a_inf<_REAL>(func, x);
}

//Calculations of other parts of the pair production rate [I]
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_TT_integrand(_REAL chi_phot, _REAL chi_ele) const
{
    if (chi_ele == zero)
        return zero;
    _REAL xx = compute_x(chi_phot, chi_ele);
    _REAL xx_3_2 = pow(xx, three/two);

    _REAL div = static_cast<_REAL>(pi)*sqrt(three);

    return compute_inner_integral(xx)-(two-chi_phot*xx_3_2)
        *k_v(two/three, (two/three)*xx_3_2)/div;
}

//Calculations of other parts of the pair production rate [II]
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_TT_function(_REAL chi_phot) const
{

    auto func = [chi_phot, this](_REAL chi_ele){
        if(chi_ele - chi_phot == zero || chi_ele == zero)
            return zero;
        else
            return compute_TT_integrand(chi_phot, chi_ele);
    };

    return quad_a_b<_REAL>(func, zero, chi_phot);
}

//Computes the cumulative pair production rate given
//chi_phot and chi_ele
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_cumulative_pair(_REAL chi_phot, _REAL chi_ele) const
{
    auto func = [this, chi_phot](_REAL chi_ele){
        return compute_TT_integrand(chi_phot, chi_ele);
    };
   _REAL num = quad_a_b<_REAL>(func, zero, chi_ele);
   return num/compute_TT_function(chi_phot) ;
}


#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
