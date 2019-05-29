#ifndef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
#define __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__

//This .hpp file contais the implementation of the
//nonlinear Breit Wheeler engine

#include<limits>

//Uses openMP to speed up the generation of the lookup table
#include <omp.h>

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

//Uses picsar vectors
#include "picsar_vector.hpp"

//Uses picsar arrays
#include "picsar_array.hpp"

//############################################### Declaration

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


        //This struct holds all the data required to re-generate
        //the breit_wheeler_engine object
        template<typename _REAL, class _RNDWRAP>
        struct breit_wheeler_innards{
                breit_wheeler_engine_ctrl<_REAL> bw_ctrl;
                _REAL lambda;
                _RNDWRAP* rng_ptr;
                size_t TTfunc_table_coords_how_many;
                _REAL* TTfunc_table_coords_ptr;
                size_t TTfunc_table_data_how_many;
                _REAL* TTfunc_table_data_ptr;
                size_t cum_distrib_table_coords_1_how_many;
                _REAL* cum_distrib_table_coords_1_ptr;
                size_t cum_distrib_table_coords_2_how_many;
                _REAL* cum_distrib_table_coords_2_ptr;
                size_t cum_distrib_table_data_how_many;
                _REAL* cum_distrib_table_data_ptr;
        };

      //Templates are used for the numerical type and for the
      //RNG wrapper
     template<typename _REAL, class _RNDWRAP>
     class breit_wheeler_engine
     {
     public:
         //A random number generator has to be passed by move.
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
         breit_wheeler_engine(const breit_wheeler_engine& other);

         //Move constructor
         breit_wheeler_engine(breit_wheeler_engine&& other);

         //Getter & setter for lambda
         _REAL get_lambda() const;
         void set_lambda(_REAL lambda);

         //get a single optical depth
         PXRMP_FORCE_INLINE
         _REAL get_optical_depth();

         //______________________GPU
         //get a single optical depth
         //same as above but conceived for GPU usage
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static _REAL internal_get_optical_depth(_REAL unf_zero_one_minus_epsi);
         //______________________

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

         //______________________GPU
         //Interp the dN_dt from table (with GPU use directly this one!)
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static _REAL
         internal_interp_dN_dt
         (_REAL energy_phot, _REAL chi_phot,
         const lookup_1d<_REAL>& ref_TTfunc_table,
         const breit_wheeler_engine_ctrl<_REAL>& ref_bw_ctrl,
         _REAL _lambda);
         //______________________

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

         //______________________GPU
         //Same as above, but with GPU use directly this one!
         //returs false if errors occur
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static bool
         internal_evolve_opt_depth_and_determine_event(
         _REAL px, _REAL py, _REAL pz,
         _REAL ex, _REAL ey, _REAL ez,
         _REAL bx, _REAL by, _REAL bz,
         _REAL dt, _REAL& opt_depth,
         bool& has_event_happend,  _REAL& event_dt,
         _REAL _lambda ,
         const lookup_1d<_REAL>& ref_TTfunc_table,
         const breit_wheeler_engine_ctrl<_REAL>& ref_bw_ctrl);
         //______________________

         //Computes the cumulative pair production rate given
         //chi_phot and chi_part
         PXRMP_FORCE_INLINE
         _REAL compute_cumulative_pair (_REAL chi_phot, _REAL chi_part) const;

         //Computes the cumulative pair production rate lookup table
         void compute_cumulative_pair_table (std::ostream* stream = nullptr);

         //This function computes the properties of the electron-positron pairs
         //generated in a BW process.
         //It is intended to be used as follows:
         //auto all = bw_engine.generate_breit_wheeler_pairs
         //   (mom[0], mom[1], mom[2],
         //    field[0], field[1], field[2],
         //    field[3], field[4], field[5], weight, sampling);
         //auto all_electrons = all[0];
         //auto all_positrons = all[1];
         //At this point all_electrons would be a vector
         //of pairs (momentum, new_weight), where new_weight
         //is simply weight/sampling.
         PXRMP_FORCE_INLINE
         picsar_array<picsar_vector<std::pair<vec3<_REAL>, _REAL>>,2>
          generate_breit_wheeler_pairs(
         _REAL px, _REAL py, _REAL pz,
         _REAL ex, _REAL ey, _REAL ez,
         _REAL bx, _REAL by, _REAL bz,
         _REAL weight, size_t sampling) ;

         //______________________GPU
         //Same as above, but with GPU use directly this one!
         //Returns false if errors occur
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static
         bool
         internal_generate_breit_wheeler_pairs(
        _REAL px, _REAL py, _REAL pz,
        _REAL ex, _REAL ey, _REAL ez,
        _REAL bx, _REAL by, _REAL bz,
        _REAL weight, size_t sampling,
        _REAL* e_px, _REAL* e_py, _REAL* e_pz,
        _REAL* p_px, _REAL* p_py, _REAL* p_pz,
        _REAL* e_weight, _REAL* p_weight,
        _REAL lambda ,
        const picsar::multi_physics::lookup_2d<_REAL>& ref_cum_distrib_table,
        const picsar::multi_physics::breit_wheeler_engine_ctrl<_REAL>& ref_bw_ctrl,
        _REAL* unf_zero_one_minus_epsi);
         //______________________

         //Write pair production table to disk
         void write_dN_dt_table(std::string filename);

         //Write cumulative_pair_table to disk
         void write_cumulative_pair_table(std::string filename);

         //Read pair production table from disk (Warning,
         //breit_wheeler_engine_ctrl is not changed in current implementation)
         void read_dN_dt_table(std::string filename);

         //Read cumulative_pair_table table from disk (Warning,
         //breit_wheeler_engine_ctrl is not changed in current implementation)
         void read_cumulative_pair_table(std::string filename);

         //Export innards
         breit_wheeler_innards<_REAL, _RNDWRAP> export_innards();

         //Constructor using innards
         breit_wheeler_engine
         (breit_wheeler_innards<_REAL, _RNDWRAP> innards);

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

        //Some handy constants
        static constexpr _REAL zero = static_cast<_REAL>(0.0);
        static constexpr _REAL one = static_cast<_REAL>(1.0);
        static constexpr _REAL two = static_cast<_REAL>(2.0);
        static constexpr _REAL three = static_cast<_REAL>(3.0);

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
breit_wheeler_engine(const breit_wheeler_engine& other):
    lambda(other.lambda), rng(other.rng), bw_ctrl(other.bw_ctrl),
    TTfunc_table(other.TTfunc_table), cum_distrib_table(other.cum_distrib_table)
    {}

//Move constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
breit_wheeler_engine(breit_wheeler_engine&& other):
    lambda(std::move(other.lambda)), rng(std::move(other.rng)),
    bw_ctrl(std::move(other.bw_ctrl)),
    TTfunc_table(std::move(other.TTfunc_table)),
    cum_distrib_table(std::move(other.cum_distrib_table))
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
    return rng.exp(one);
}

//______________________GPU
//get a single optical depth
//same as above but conceived for GPU usage
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
internal_get_optical_depth(_REAL unf_zero_one_minus_epsi)
{
    return -log(one - unf_zero_one_minus_epsi);
}
//______________________

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
    //Prepare the TT_coords vector
    picsar_vector<_REAL> TT_coords = generate_log_spaced_vec(
     bw_ctrl.chi_phot_tdndt_min, bw_ctrl.chi_phot_tdndt_max,
    bw_ctrl.chi_phot_tdndt_how_many);

    msg("Computing table for dNdt...\n", stream);
    //Do the hard work
    picsar_vector<_REAL> TT_vals{TT_coords.size()};
    std::transform(TT_coords.begin(), TT_coords.end(),
        TT_vals.begin(),
        [this](_REAL chi){return compute_TT_function(chi);});
    msg("...done!\n", stream);

    //The table will store the logarithms
    auto logfun = [](_REAL val){return log(val);};
    std::transform(TT_coords.begin(), TT_coords.end(), TT_coords.begin(), logfun);
    std::transform(TT_vals.begin(), TT_vals.end(), TT_vals.begin(), logfun);

    TTfunc_table = lookup_1d<_REAL>{TT_coords, TT_vals};
}

//Interp the dN_dt from table
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
interp_dN_dt(_REAL energy_phot, _REAL chi_phot) const
{
    return internal_interp_dN_dt
        (energy_phot, chi_phot, TTfunc_table, bw_ctrl, lambda);
}


//Interp the dN_dt from table (with GPU use directly this one!)
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
internal_interp_dN_dt
(_REAL energy_phot, _REAL chi_phot,
const picsar::multi_physics::lookup_1d<_REAL>& ref_TTfunc_table,
const picsar::multi_physics::breit_wheeler_engine_ctrl<_REAL>& ref_bw_ctrl,
_REAL _lambda)
{
#ifdef PXRMP_WITH_SI_UNITS
    _lambda = static_cast<_REAL>(1.0);
#endif


    if(energy_phot == zero || chi_phot == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__pair_prod_coeff)*
        _lambda*(one/( chi_phot * energy_phot));

        _REAL TT = zero;
    //Use approximate analytical expression if chi < chi_phot_tdndt_min
    //or chi > chi_phot_tdndt_min

    _REAL a = static_cast<_REAL>(__erber_Tfunc_asynt_a);
    _REAL b = static_cast<_REAL>(__erber_Tfunc_asynt_b);
    if(chi_phot <= ref_bw_ctrl.chi_phot_tdndt_min){
        //Not suitable for GPU! Replaced with asymptotic expansion.
        //TT = a*chi_phot*pow(k_v(one/three, b/chi_phot),two);
        TT = (static_cast<_REAL>(pi)*a/(two*b))*chi_phot*chi_phot*
            exp(-two*b/chi_phot);
    }
    else if(chi_phot >= ref_bw_ctrl.chi_phot_tdndt_max){
        //Not suitable for GPU! Replaced with asymptotic expansion.
        //TT = a*chi_phot*pow(k_v(one/three, b/chi_phot),two);
        _REAL ccoeff = static_cast<_REAL>(tgamma(one/three)/two);
        TT = a*chi_phot*ccoeff*pow(chi_phot*two/b , two/three);

    }
    //otherwise use lookup tables
    else{
            TT =  exp(ref_TTfunc_table.interp_linear_equispaced(log(chi_phot)));
    }
    //**end

    _REAL dndt = coeff * TT;

    return dndt;
}


//This function evolves the optical depth for a particle and
//checks if it goes to zero. If it doesn't the output is false,0.
//On the contrary, if it goes to zero the output is true, dt_em,
//where dt_em (<= dt) is the time at which the event occurs.
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
    bool has_event_happend = false;
    _REAL event_dt = zero;

    if(!internal_evolve_opt_depth_and_determine_event(px, py, pz, ex, ey, ez,
        bx, by, bz, dt, opt_depth, has_event_happend, event_dt, lambda,
        TTfunc_table, bw_ctrl)){

        err("dndt lookup table not initialized!\n");

        _REAL energy = norm<_REAL>(vec3<_REAL> {px, py, pz})*__c;
        _REAL chi = chi_photon(px, py, pz, ex, ey, ez, bx, by, bz, lambda);
        _REAL dndt = compute_dN_dt(energy, chi);
        opt_depth -= dndt*dt;
        if(opt_depth < zero){
            //Calculates the time at which pair prod takes place
            event_dt = opt_depth/dndt + dt;
            has_event_happend = true;
        }
    }

    return std::make_pair(has_event_happend, event_dt);
}

//Same as above but tailored for direct use with GPU
//returns false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
bool
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
internal_evolve_opt_depth_and_determine_event(
_REAL px, _REAL py, _REAL pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dt, _REAL& opt_depth,
bool& has_event_happend, _REAL& event_dt,
_REAL _lambda,
const picsar::multi_physics::lookup_1d<_REAL>& ref_TTfunc_table,
const picsar::multi_physics::breit_wheeler_engine_ctrl<_REAL>& ref_bw_ctrl)
{

#ifdef PXRMP_WITH_SI_UNITS
    _lambda = static_cast<_REAL>(1.0);
#endif

    _REAL energy = norm<_REAL>(vec3<_REAL> {px, py, pz})*__c;
    _REAL chi = chi_photon(px, py, pz, ex, ey, ez, bx, by, bz, _lambda);

    has_event_happend = false;
    event_dt = zero;

    //Do NOT evolve opt_depth if the chi parameter is less then threshold
    //or if the photon energy is not high enough to generate a pair
    if(chi <= ref_bw_ctrl.chi_phot_min ||
       energy < two*static_cast<_REAL>(__emass*__c*__c))
           return true;

    //**Compute dndt
    _REAL dndt;
    //Uses table if available
    if(ref_TTfunc_table.is_init()){
            dndt =
            internal_interp_dN_dt(energy, chi, ref_TTfunc_table, ref_bw_ctrl, _lambda);
    }
    //If not..
    else{
        return false;
    }

    opt_depth -= dndt*dt;

    if(opt_depth < zero){
            //Calculates the time at which pair prod takes place
        event_dt = opt_depth/dndt + dt;
        has_event_happend = true;
    }

    return true;
}


//Computes the cumulative pair production rate given
//chi_phot and chi_part
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_cumulative_pair(_REAL chi_phot, _REAL chi_part) const
{
    auto func = [this, chi_phot](_REAL chi_part){
        return compute_TT_integrand(chi_phot, chi_part);
    };
   _REAL num = quad_a_b<_REAL>(func, zero, chi_part);
   return num/compute_TT_function(chi_phot) ;
}



template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
compute_cumulative_pair_table (std::ostream* stream)
{
    //Prepare the chi_coords vector
    picsar_vector<_REAL> chi_coords = generate_log_spaced_vec(
     bw_ctrl.chi_phot_tpair_min, bw_ctrl.chi_phot_tpair_max,
    bw_ctrl.chi_phot_tpair_how_many);

    picsar_vector<_REAL> frac_coords = generate_lin_spaced_vec(zero,one/two,
    bw_ctrl.chi_frac_tpair_how_many);

    picsar_vector<_REAL> pair_vals{chi_coords.size()*frac_coords.size()};

    msg("Computing table for pair production...\n", stream);



    #pragma omp parallel for
    for(size_t ii = 0; ii < chi_coords.size(); ii++){
        auto chi_phot = chi_coords[ii];
        size_t cc = ii*frac_coords.size();
        pair_vals[cc++] = zero;
        for(size_t jj = 1; jj < frac_coords.size() - 1; jj++){
            pair_vals[cc++] = compute_cumulative_pair(
                chi_phot, chi_phot*frac_coords[jj]);
        }
        pair_vals[cc++] = (one/two); //The function is symmetric
    }
    msg("...done!\n", stream);

    //The table will store the logarithms
    auto logfun = [](_REAL val){return log(val);};
    std::transform(chi_coords.begin(), chi_coords.end(), chi_coords.begin(), logfun);
    std::transform(pair_vals.begin(), pair_vals.end(), pair_vals.begin(), logfun);

    cum_distrib_table = lookup_2d<_REAL>{
        picsar_array<picsar_vector<_REAL>,2>{chi_coords, frac_coords}, pair_vals};

}

//This function computes the properties of the electron-positron pairs
//generated in a BW process.
//It is intended to be used as follows:
//auto all = bw_engine.generate_breit_wheeler_pairs
//   (mom[0], mom[1], mom[2],
//    field[0], field[1], field[2],
//    field[3], field[4], field[5], weight, sampling);
//auto all_electrons = all[0];
//auto all_positrons = all[1];
//At this point all_electrons would be a vector
//of pairs (momentum, new_weight), where new_weight
//is simply weight/sampling.
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_array<
picsar::multi_physics::picsar_vector<
std::pair<picsar::multi_physics::vec3<_REAL>, _REAL>>,2>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
generate_breit_wheeler_pairs(
_REAL px, _REAL py, _REAL pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL weight, size_t sampling)
{
    picsar_vector<_REAL> e_px{sampling};
    picsar_vector<_REAL> e_py{sampling};
    picsar_vector<_REAL> e_pz{sampling};
    picsar_vector<_REAL> p_px{sampling};
    picsar_vector<_REAL> p_py{sampling};
    picsar_vector<_REAL> p_pz{sampling};
    picsar_vector<_REAL> e_weight{sampling};
    picsar_vector<_REAL> p_weight{sampling};

    picsar_vector<_REAL> unf_zero_one_minus_epsi{sampling};
    for(auto& el: unf_zero_one_minus_epsi)
        el = rng.unf(zero, one);

    //Call to the static GPU-friendly function
    internal_generate_breit_wheeler_pairs(
    px, py, pz, ex, ey, ez, bx, by, bz, weight, sampling,
    e_px.data(), e_py.data(), e_pz.data(),
    p_px.data(), p_py.data(), p_pz.data(),
    e_weight.data(), p_weight.data(),
    lambda, cum_distrib_table, bw_ctrl,
    unf_zero_one_minus_epsi.data());

    picsar_vector<std::pair<vec3<_REAL>, _REAL>> electrons(sampling);
    picsar_vector<std::pair<vec3<_REAL>, _REAL>> positrons(sampling);

    for(size_t s = 0; s < sampling; s++){
        electrons[s] =
        std::make_pair(vec3<_REAL>{e_px[s], e_py[s], e_pz[s]}, e_weight[s]);
        positrons[s] =
        std::make_pair(vec3<_REAL>{p_px[s], p_py[s], p_pz[s]}, p_weight[s]);
    }

    return {electrons, positrons};
}

//______________________GPU
//Same as above, but with GPU use directly this one!
//Returns false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
bool
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
internal_generate_breit_wheeler_pairs(
_REAL px, _REAL py, _REAL pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL weight, size_t sampling,
_REAL* e_px, _REAL* e_py, _REAL* e_pz,
_REAL* p_px, _REAL* p_py, _REAL* p_pz,
_REAL* e_weight, _REAL* p_weight,
_REAL _lambda ,
const picsar::multi_physics::lookup_2d<_REAL>& ref_cum_distrib_table,
const picsar::multi_physics::breit_wheeler_engine_ctrl<_REAL>& ref_bw_ctrl,
_REAL* unf_zero_one_minus_epsi)
{
    _REAL chi_phot = chi_photon(px, py, pz, ex, ey, ez, bx, by, bz, _lambda);

    _REAL new_weight = weight/sampling;

    _REAL me_c = static_cast<_REAL>(__emass*__c);

    vec3<_REAL> p_phot{px, py, pz};
    _REAL norm_phot = norm(p_phot);
    vec3<_REAL> n_phot = p_phot/norm_phot;
    _REAL gamma_phot = norm_phot/me_c;

    const size_t how_many_frac = ref_cum_distrib_table.ref_coords()[1].size();

    _REAL tab_chi_phot = chi_phot;
    if(chi_phot < ref_bw_ctrl.chi_phot_tpair_min)
        tab_chi_phot = ref_bw_ctrl.chi_phot_tpair_min;
    else if(chi_phot > ref_bw_ctrl.chi_phot_tpair_max)
        tab_chi_phot = ref_bw_ctrl.chi_phot_tpair_max;

    for(size_t s = 0; s < sampling; s++){
        _REAL prob = unf_zero_one_minus_epsi[s];
        bool invert = false;
        if(prob >= one/two){
            prob -= one/two;
            invert = true;
        }

        size_t first = 0;
        size_t it;
        _REAL val;
        size_t count = how_many_frac;
        while(count > 0){
            it = first;
            size_t step = count/2;
            it += step;
            val =
                exp(ref_cum_distrib_table.interp_linear_first_equispaced
                    (log(tab_chi_phot), it));
            if(!(prob < val)){
                first = ++it;
                count -= step+1;
            }
            else{
                count = step;
            }
        }

        size_t upper = first;
        size_t lower = upper-1;

        const _REAL upper_frac = ref_cum_distrib_table.ref_coords()[1][upper];
        const _REAL lower_frac = ref_cum_distrib_table.ref_coords()[1][lower];
        _REAL upper_prob =
            exp(ref_cum_distrib_table.interp_linear_first_equispaced
            (log(tab_chi_phot), upper));
        _REAL lower_prob =
            exp(ref_cum_distrib_table.interp_linear_first_equispaced
            (log(tab_chi_phot), lower));
        _REAL chi_ele_frac = lower_frac +
            (prob-lower_prob)*(upper_frac-lower_frac)/(upper_prob-lower_prob);

        _REAL chi_ele = chi_ele_frac*chi_phot;
        _REAL chi_pos = chi_phot - chi_ele;

        if(invert){
            _REAL tt = chi_ele;
            chi_ele = chi_pos;
            chi_pos = tt;
        }

        _REAL coeff =  (gamma_phot - two)/chi_phot;

        _REAL cc_ele = sqrt((one+chi_ele*coeff)*(one+chi_ele*coeff)-one)*me_c;
        _REAL cc_pos = sqrt((one+chi_pos*coeff)*(one+chi_pos*coeff)-one)*me_c;

        vec3<_REAL> p_ele = cc_ele*n_phot;
        vec3<_REAL> p_pos = cc_pos*n_phot;


        e_px[s] =  p_ele[0];
        e_py[s] =  p_ele[1];
        e_pz[s] =  p_ele[2];

        p_px[s] =  p_pos[0];
        p_py[s] =  p_pos[1];
        p_pz[s] =  p_pos[2];

        e_weight[s] = new_weight;
        p_weight[s] = new_weight;
    }

    return true;
}

//Write pair production table to disk
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
write_dN_dt_table(std::string filename)
{
    std::ofstream of;
    of.open(filename, std::ios::binary);
    TTfunc_table.write_on_stream_bin(of);
    of.close();
}

//Write cumulative_pair_table to disk
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
write_cumulative_pair_table(std::string filename)
{
    std::ofstream of;
    of.open(filename, std::ios::binary);
    cum_distrib_table.write_on_stream_bin(of);
    of.close();
}

//Read pair production table from disk (Warning,
//breit_wheeler_engine_ctrl is not changed in current implementation)
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
read_dN_dt_table(std::string filename)
{
    std::ifstream iif;
    iif.open(filename, std::ios::binary);
    TTfunc_table.read_from_stream_bin(iif);
    iif.close();
}

//Read cumulative_pair_table table from disk (Warning,
//breit_wheeler_engine_ctrl is not changed in current implementation)
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::
read_cumulative_pair_table(std::string filename)
{
    std::ifstream iif;
    iif.open(filename, std::ios::binary);
    cum_distrib_table.read_from_stream_bin(iif);
    iif.close();
}

//Export innards
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_innards<_REAL, _RNDWRAP>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::export_innards()
{
    breit_wheeler_innards<_REAL, _RNDWRAP> innards;

    //Filling innards
    innards.bw_ctrl                         = bw_ctrl;
    innards.lambda                          = lambda;
    innards.rng_ptr                         = &rng;
    innards.TTfunc_table_coords_how_many    = TTfunc_table.ref_coords().size();
    innards.TTfunc_table_coords_ptr         = TTfunc_table.ref_coords().data();
    innards.TTfunc_table_data_how_many      = TTfunc_table.ref_data().size();
    innards.TTfunc_table_data_ptr           = TTfunc_table.ref_data().data();
    innards.cum_distrib_table_coords_1_how_many
        = cum_distrib_table.ref_coords()[0].size();
    innards.cum_distrib_table_coords_1_ptr
        = cum_distrib_table.ref_coords()[0].data();
    innards.cum_distrib_table_coords_2_how_many
        = cum_distrib_table.ref_coords()[1].size();
    innards.cum_distrib_table_coords_2_ptr
            = cum_distrib_table.ref_coords()[1].data();
    innards.cum_distrib_table_data_how_many
        = cum_distrib_table.ref_data().size();
    innards.cum_distrib_table_data_ptr
        = cum_distrib_table.ref_data().data();

    return innards;
}

//Constructor using innards
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::breit_wheeler_engine<_REAL, _RNDWRAP>::breit_wheeler_engine
(picsar::multi_physics::breit_wheeler_innards<_REAL, _RNDWRAP> innards):
    bw_ctrl{innards.bw_ctrl},
    lambda{innards.lambda},
    TTfunc_table{innards.TTfunc_table_coords_how_many,
                 innards.TTfunc_table_coords_ptr,
                 innards.TTfunc_table_data_ptr},
    cum_distrib_table{innards.cum_distrib_table_coords_1_how_many,
                      innards.cum_distrib_table_coords_1_ptr,
                      innards.cum_distrib_table_coords_2_how_many,
                      innards.cum_distrib_table_coords_2_ptr,
                      innards.cum_distrib_table_data_ptr},
    rng{*innards.rng_ptr}
{}

//___________________PRIVATE FUNCTIONS_____________________________________

//Function to compute X (Warning: it doesn't chek if chi_ele != 0 or
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
    auto func = [=](_REAL s){
        if (s >= static_cast<_REAL>(__breit_wheeler_special_func_big_arg))
            return zero;
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
    if (chi_phot == zero || chi_ele >= chi_phot || chi_ele == zero)
        return zero;
    _REAL xx = compute_x(chi_phot, chi_ele);
    _REAL xx_3_2 = pow(xx, three/two);

    _REAL div = static_cast<_REAL>(pi)*sqrt(three);

    return (compute_inner_integral(xx)-(two-chi_phot*xx_3_2)
        *k_v(two/three, (two/three)*xx_3_2))/div;
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

#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
