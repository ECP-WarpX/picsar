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

//Uses picsar vectors
#include "picsar_vector.hpp"

//Uses picsar arrays
#include "picsar_array.hpp"

//Uses utilities
#include "utilities.hpp"

//############################################### Declaration

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
          _REAL chi_frac_tem_min =
            static_cast<_REAL>(__quantum_synchrotron_chi_frac_tem_min_frac);
          size_t chi_frac_tem_how_many =
            __quantum_synchrotron_chi_frac_tem_how_many;

         //Minimum chi_phot to consider in integrations:
         _REAL chi_min = chi_frac_tem_min*chi_part_tem_min;
      };

      //This struct holds all the data required to re-generate
      //the quantum_synchrotron_engine object
      template<typename _REAL, class _RNDWRAP>
      struct quantum_synchrotron_innards{
              quantum_synchrotron_engine_ctrl<_REAL> qs_ctrl;
              _REAL lambda;
              _RNDWRAP* rng_ptr;
              size_t KKfunc_table_coords_how_many;
              _REAL* KKfunc_table_coords_ptr;
              size_t KKfunc_table_data_how_many;
              _REAL* KKfunc_table_data_ptr;
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
         quantum_synchrotron_engine_ctrl<_REAL> qs_ctrl =
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

         //______________________GPU
         //get a single optical depth
         //same as above but conceived for GPU usage
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static _REAL internal_get_optical_depth(_REAL unf_zero_one_minus_epsi);


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

         //______________________GPU
         //Interp the dN_dt from table (with GPU use directly this one!)
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static _REAL
         internal_interp_dN_dt
         (_REAL energy_part, _REAL chi_part,
         const lookup_1d<_REAL>& ref_KKfunc_table,
         const quantum_synchrotron_engine_ctrl<_REAL>& ref_qs_ctrl,
         _REAL _lambda);
         //______________________


         //This function evolves the optical depth for a particle and
         //checks if it goes to zero. If it doesn't the output is false,0.
         //On the contrary, if it goes to zero the output is true, dt_phot,
         //where dt_phot (<= dt) is the time at which the event occurs.
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
         const lookup_1d<_REAL>& ref_KKfunc_table,
         const quantum_synchrotron_engine_ctrl<_REAL>& ref_qs_ctrl);
         //______________________


         //Computes the cumulative photon emission rate given
         //chi_phot and chi_part
         PXRMP_FORCE_INLINE
         _REAL compute_cumulative_phot_em(_REAL chi_phot, _REAL chi_part) const;

         //Computes the cumulative photon emission rate lookup table
         void compute_cumulative_phot_em_table (std::ostream* stream = nullptr);

         //This function computes the properties of the photons
         //generated in a quantum synchrotron-like process.
         //It is intended to be used as follows:
         //auto all_photons = qs_engine.generate_photons
         //   (mom[0], mom[1], mom[2],
         //    field[0], field[1], field[2],
         //    field[3], field[4], field[5], weight, sampling);
         //At this point all_photons would be a vector
         //of pairs (momentum, new_weight), where new_weight
         //is simply weight/sampling.
         //The function updates also the momentum of the particle
         PXRMP_FORCE_INLINE
         picsar_vector<std::pair<vec3<_REAL>, _REAL>>
         generate_photons_and_update_momentum(
         _REAL& px, _REAL& py, _REAL& pz,
         _REAL ex, _REAL ey, _REAL ez,
         _REAL bx, _REAL by, _REAL bz,
         _REAL weight, size_t sampling);

         //______________________GPU
         //Same as above, but with GPU use directly this one!
         //Returns false if errors occur
         PXRMP_GPU
         PXRMP_FORCE_INLINE
         static
         bool
         internal_generate_photons_and_update_momentum(
         _REAL& px, _REAL& py, _REAL& pz,
         _REAL ex, _REAL ey, _REAL ez,
         _REAL bx, _REAL by, _REAL bz,
         _REAL weight, size_t sampling,
         _REAL* g_px, _REAL* g_py, _REAL* g_pz,
         _REAL* g_weight,
         _REAL lambda ,
         const picsar::multi_physics::lookup_2d<_REAL>& ref_cum_distrib_table,
         const picsar::multi_physics::quantum_synchrotron_engine_ctrl<_REAL>& ref_qs_ctrl,
         _REAL* unf_zero_one_minus_epsi
         );
         //______________________

         //Write total cross section table to disk
         void write_dN_dt_table(std::string filename);

         //Write cumulative photon emission rate to disk
         void write_cumulative_phot_em_table(std::string filename);

         //Read total cross section table from disk (Warning,
         //qs_ctrl is not changed in current implementation)
         void read_dN_dt_table(std::string filename);

         //Read cumulative photon emission rate to disk (Warning,
         //qs_ctrl is not changed in current implementation)
         void read_cumulative_phot_em_table(std::string filename);

         //Export innards
         quantum_synchrotron_innards<_REAL, _RNDWRAP> export_innards();

         //Constructor using innards
         quantum_synchrotron_engine
         (quantum_synchrotron_innards<_REAL, _RNDWRAP> innards);



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

         //Some handy constants
         static constexpr _REAL zero = static_cast<_REAL>(0.0);
         static constexpr _REAL one = static_cast<_REAL>(1.0);
         static constexpr _REAL two = static_cast<_REAL>(2.0);
         static constexpr _REAL three = static_cast<_REAL>(3.0);
         static constexpr _REAL four = static_cast<_REAL>(4.0);
         static constexpr _REAL five = static_cast<_REAL>(5.0);

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
    KKfunc_table(other.KKfunc_table), cum_distrib_table(other.cum_distrib_table)
    {}

//Move constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
quantum_synchrotron_engine(quantum_synchrotron_engine&& other):
    lambda(std::move(other.lambda)), rng(std::move(other.rng)),
    qs_ctrl(std::move(other.qs_ctrl)),
    KKfunc_table(std::move(other.KKfunc_table)),
    cum_distrib_table(std::move(other.cum_distrib_table))
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

//______________________GPU
//get a single optical depth
//same as above but conceived for GPU usage
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
internal_get_optical_depth(_REAL unf_zero_one_minus_epsi)
{
    return -log(one - unf_zero_one_minus_epsi);
}
//______________________


//Calculates the photon emission rate (Warning: no lookup tables)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_dN_dt(_REAL energy_part, _REAL chi_part) const
{
    if(energy_part == zero || chi_part == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__quantum_synchrotron_rate_coeff)*
    lambda*(one/energy_part);

    return coeff*compute_KK_function(chi_part);
}

//Computes the lookup_table needed for dN/dt
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
void picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_dN_dt_lookup_table(std::ostream* stream)
{
    //Prepare the KK_coords vector
    picsar_vector<_REAL> KK_coords = generate_log_spaced_vec(
     qs_ctrl.chi_part_tdndt_min, qs_ctrl.chi_part_tdndt_max,
    qs_ctrl.chi_part_tdndt_how_many);

    msg("Computing table for dNdt...\n", stream);
    //Do the hard work
    picsar_vector<_REAL> KK_vals{KK_coords.size()};
    std::transform(KK_coords.begin(), KK_coords.end(),
        KK_vals.begin(),
        [this](_REAL chi){return compute_KK_function(chi);});
    msg("...done!\n", stream);

    //The table will store the logarithms
    auto logfun = [](_REAL val){return log(val);};
    std::transform(KK_coords.begin(), KK_coords.end(), KK_coords.begin(), logfun);
    std::transform(KK_vals.begin(), KK_vals.end(), KK_vals.begin(), logfun);

    KKfunc_table = lookup_1d<_REAL>{KK_coords, KK_vals};
}

//Interp the dN_dt from table
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
interp_dN_dt(_REAL energy_part, _REAL chi_part) const
{
    return internal_interp_dN_dt
        (energy_part, chi_part, KKfunc_table, qs_ctrl, lambda);
}

//Interp the dN_dt from table (with GPU use directly this one!)
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
internal_interp_dN_dt
(_REAL energy_part, _REAL chi_part,
const lookup_1d<_REAL>& ref_KKfunc_table,
const quantum_synchrotron_engine_ctrl<_REAL>& ref_qs_ctrl,
_REAL _lambda)
{
    if(energy_part == zero || chi_part == zero)
        return zero;

    _REAL coeff = static_cast<_REAL>(__quantum_synchrotron_rate_coeff)*
        _lambda*one/energy_part;

    _REAL KK = zero;

    //If chi is out of table, use the first (or the last) value
    //in the table
    if(chi_part <= ref_qs_ctrl.chi_part_tdndt_min){
        chi_part = ref_qs_ctrl.chi_part_tdndt_min;
    }
    else if(chi_part >= ref_qs_ctrl.chi_part_tdndt_max){
        chi_part = ref_qs_ctrl.chi_part_tdndt_max ;
    }

    KK =  exp(ref_KKfunc_table.interp_linear_equispaced(log(chi_part)));

    _REAL dndt = coeff * KK;

    return dndt;
}

//This function evolves the optical depth for a particle and
//checks if it goes to zero. If it doesn't the output is false,0.
//On the contrary, if it goes to zero the output is true, dt_phot,
//where dt_phot (<= dt) is the time at which the event occurs.
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
std::pair<bool, _REAL>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
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
        KKfunc_table, qs_ctrl)){

        err("dndt lookup table not initialized!\n");

        _REAL energy = norm<_REAL>(vec3<_REAL> {px, py, pz})*__c;
        _REAL chi = chi_lepton(px, py, pz, ex, ey, ez, bx, by, bz, lambda);
        _REAL dndt = compute_dN_dt(energy, chi);
        opt_depth -= dndt*dt;
        if(opt_depth < zero){
            //Calculates the time at which photon emission takes place
            event_dt = opt_depth/dndt + dt;
            has_event_happend = true;
        }
    }

    return std::make_pair(has_event_happend, event_dt);
}


//______________________GPU
//Same as above, but with GPU use directly this one!
//returs false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
bool
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
internal_evolve_opt_depth_and_determine_event(
_REAL px, _REAL py, _REAL pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dt, _REAL& opt_depth,
bool& has_event_happend,  _REAL& event_dt,
_REAL _lambda ,
const lookup_1d<_REAL>& ref_KKfunc_table,
const quantum_synchrotron_engine_ctrl<_REAL>& ref_qs_ctrl)
{
    _REAL energy = norm<_REAL>(vec3<_REAL> {px, py, pz})*__c;
    _REAL chi = chi_lepton(px, py, pz, ex, ey, ez, bx, by, bz, _lambda);

    has_event_happend = false;
    event_dt = zero;

    //Do NOT evolve opt_depth if the chi parameter is less then threshold
    if(chi <= ref_qs_ctrl.chi_part_min)
        return true;

    //**Compute dndt
    _REAL dndt;
    //Uses table if available
    if(ref_KKfunc_table.is_init()){
            dndt =
                internal_interp_dN_dt
                (energy, chi, ref_KKfunc_table, ref_qs_ctrl, _lambda);
    }
    //If not..
    else{
            return false;
    }

    opt_depth -= dndt*dt;

    if(opt_depth < zero){
        //Calculates the time at which photon emission takes place
        event_dt = opt_depth/dndt + dt;
        has_event_happend = true;
    }

    return true;
}
//______________________


//Computes the cumulative photon emission rate given
//chi_phot and chi_part
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_cumulative_phot_em(_REAL chi_phot, _REAL chi_part) const
{
    auto func = [chi_part, this](_REAL __chi_phot){

        if(chi_part - __chi_phot == zero || __chi_phot == zero)
            return zero;
        else{
            return compute_KK_integrand(chi_part, __chi_phot);
        }

    };
    _REAL num = quad_a_b<_REAL>(func, zero, chi_phot/two) + quad_a_b<_REAL>(func, chi_phot/two, chi_phot);
    _REAL kk = compute_KK_function(chi_part);
    if(num >= kk)
        return one;
    return num/kk;
}

//Computes the cumulative photon emission rate lookup table
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
void
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_cumulative_phot_em_table (std::ostream* stream)
{
    //Prepare the chi_coords vector
    picsar_vector<_REAL> chi_coords = generate_log_spaced_vec(
     qs_ctrl.chi_part_tem_min, qs_ctrl.chi_part_tem_max,
    qs_ctrl.chi_part_tem_how_many);

    picsar_vector<_REAL> frac_coords = generate_log_spaced_vec(qs_ctrl.chi_frac_tem_min, one,
    qs_ctrl.chi_frac_tem_how_many);

    picsar_vector<_REAL> phot_vals{chi_coords.size()*frac_coords.size()};

    msg("Computing table for photon emission...\n", stream);

    size_t cc = 0;
    for(auto chi_part: chi_coords){
        msg("chi_part: " + std::to_string(chi_part) + " \n", stream);
        for(size_t i = 0; i < frac_coords.size() - 1; i++){
            msg("   log10 frac: " + std::to_string(log10(frac_coords[i])) + " --> ", stream);
            phot_vals[cc++] = compute_cumulative_phot_em(
                chi_part*frac_coords[i], chi_part);
            msg(std::to_string(phot_vals[cc-1]*100) + " % \n", stream);
        }
        phot_vals[cc++] = one;
    }
    msg("...done!\n", stream);

    //The table will store the logarithms
    auto logfun = [](_REAL val){return log(val);};
    std::transform(chi_coords.begin(), chi_coords.end(), chi_coords.begin(), logfun);
    std::transform(frac_coords.begin(), frac_coords.end(), frac_coords.begin(), logfun);
    std::transform(phot_vals.begin(), phot_vals.end(), phot_vals.begin(), logfun);

    cum_distrib_table = lookup_2d<_REAL>{
        picsar_array<picsar_vector<_REAL>,2>{chi_coords, frac_coords}, phot_vals};
}


//This function computes the properties of the photons
//generated in a quantum synchrotron-like process.
//It is intended to be used as follows:
//auto all_photons = qs_engine.generate_photons
//   (mom[0], mom[1], mom[2],
//    field[0], field[1], field[2],
//    field[3], field[4], field[5], weight, sampling);
//At this point all_photons would be a vector
//of pairs (momentum, new_weight), where new_weight
//is simply weight/sampling.
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
picsar::multi_physics::
    picsar_vector<std::pair<picsar::multi_physics::vec3<_REAL>, _REAL>>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
 generate_photons_and_update_momentum(
_REAL& px, _REAL& py, _REAL& pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL weight, size_t sampling)
{
    picsar_vector<_REAL> g_px{sampling};
    picsar_vector<_REAL> g_py{sampling};
    picsar_vector<_REAL> g_pz{sampling};
    picsar_vector<_REAL> g_weight{sampling};

    picsar_vector<_REAL> unf_zero_one_minus_epsi{sampling};
    for(auto& el: unf_zero_one_minus_epsi)
        el = rng.unf(zero, one);

    //Call to the static GPU-friendly function
    internal_generate_photons_and_update_momentum(
        px, py, pz, ex, ey, ez, bx, by, bz, weight, sampling,
        g_px.data(), g_py.data(), g_pz.data(),
        g_weight.data(),
        lambda, cum_distrib_table, qs_ctrl,
        unf_zero_one_minus_epsi.data());

    picsar_vector<std::pair<vec3<_REAL>, _REAL>> photons(sampling);

    for(size_t s = 0; s < sampling; s++){
        photons[s] =
        std::make_pair(vec3<_REAL>{g_px[s], g_py[s], g_pz[s]}, g_weight[s]);
    }

    return photons;
}

//Same as above but tailored for direct use with GPU
//returns false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
bool
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
internal_generate_photons_and_update_momentum(
_REAL& px, _REAL& py, _REAL& pz,
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL weight, size_t sampling,
_REAL* g_px, _REAL* g_py, _REAL* g_pz,
_REAL* g_weight,
_REAL _lambda ,
const picsar::multi_physics::lookup_2d<_REAL>& ref_cum_distrib_table,
const picsar::multi_physics::quantum_synchrotron_engine_ctrl<_REAL>& ref_qs_ctrl,
_REAL* unf_zero_one_minus_epsi)
{

    _REAL chi_part = chi_lepton(px, py, pz, ex, ey, ez, bx, by, bz, _lambda);

    _REAL new_weight = weight/sampling;

    _REAL me_c = static_cast<_REAL>(__emass*__c);

    vec3<_REAL> p_part{px, py, pz};
    _REAL p2 = norm2(p_part);
    _REAL norm_part = sqrt(p2);
    vec3<_REAL> n_part = p_part/norm_part;
    _REAL gamma_part = sqrt(one + p2/me_c/me_c);


    const size_t how_many_frac = ref_cum_distrib_table.ref_coords()[1].size();

    _REAL tab_chi_part = chi_part;
    if(chi_part < ref_qs_ctrl.chi_part_tem_min)
        tab_chi_part = ref_qs_ctrl.chi_part_tem_min;
    else if(chi_part > ref_qs_ctrl.chi_part_tem_max)
        tab_chi_part = ref_qs_ctrl.chi_part_tem_max;

    _REAL one_over_sampling = one/sampling;

    _REAL min = exp(ref_cum_distrib_table.interp_linear_first_equispaced
            (log(tab_chi_part), 0));

    for(size_t s = 0; s < sampling; s++){
        _REAL prob = unf_zero_one_minus_epsi[s];

        _REAL chi_phot_frac;

        if(prob < min){
            chi_phot_frac = exp(ref_cum_distrib_table.ref_coords()[1][0]);
        }
        else{
            _REAL log_prob = log(prob);

            size_t first = 0;
            size_t it;
            size_t count = how_many_frac;
            while(count > 0){
                it = first;
                size_t step = count/2;
                it += step;

                _REAL val = (ref_cum_distrib_table.interp_linear_first_equispaced
                        (log(tab_chi_part), it));

                if(!(log_prob < val)){
                    first = ++it;
                    count -= step+1;
                }
                else{
                    count = step;
                }
            }

            size_t upper = first;
            size_t lower = upper-1;


            /*
            const _REAL upper_frac =
                exp(ref_cum_distrib_table.ref_coords()[1][upper]);
            const _REAL lower_frac =
                exp(ref_cum_distrib_table.ref_coords()[1][lower]);
            _REAL upper_prob =
                exp(ref_cum_distrib_table.interp_linear_first_equispaced
                (log(tab_chi_part), upper));
            _REAL lower_prob =
                exp(ref_cum_distrib_table.interp_linear_first_equispaced
                (log(tab_chi_part), lower));
            chi_phot_frac = lower_frac +
                (prob-lower_prob)*(upper_frac-lower_frac)
                /(upper_prob-lower_prob);
                */

            const _REAL upper_log_frac =
                ref_cum_distrib_table.ref_coords()[1][upper];
            const _REAL lower_log_frac =
                ref_cum_distrib_table.ref_coords()[1][lower];
            _REAL upper_log_prob =
                ref_cum_distrib_table.interp_linear_first_equispaced
                (log(tab_chi_part), upper);
            _REAL lower_log_prob =
                ref_cum_distrib_table.interp_linear_first_equispaced
                (log(tab_chi_part), lower);
            _REAL chi_phot_log_frac = lower_log_frac +
                (log_prob-lower_log_prob)*(upper_log_frac-lower_log_frac)
                /(upper_log_prob-lower_log_prob);
            chi_phot_frac = exp(chi_phot_log_frac);

        }


        _REAL chi_phot = chi_part*chi_phot_frac;

        _REAL gamma_phot = chi_phot/chi_part*(gamma_part-one);

        vec3<_REAL>  p_phot = gamma_phot * n_part * me_c;

        g_px[s] = p_phot[0];
        g_py[s] = p_phot[1];
        g_pz[s] = p_phot[2];

        g_weight[s] = new_weight;

        //Update particle momentum
        px -= p_phot[0]*one_over_sampling;
        py -= p_phot[1]*one_over_sampling;
        pz -= p_phot[2]*one_over_sampling;

    }

    return true;
}



//Write total cross section table to disk
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
write_dN_dt_table(std::string filename)
{
    std::ofstream of;
    of.open(filename, std::ios::binary);
    KKfunc_table.write_on_stream_bin(of);
    of.close();
}

//Write cumulative photon emission rate to disk
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
write_cumulative_phot_em_table(std::string filename)
{
    std::ofstream of;
    of.open(filename, std::ios::binary);
    cum_distrib_table.write_on_stream_bin(of);
    of.close();
}

//Read total cross section table from disk (Warning,
//qs_ctrl is not changed in current implementation)
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
read_dN_dt_table(std::string filename)
{
    std::ifstream iif;
    iif.open(filename, std::ios::binary);
    KKfunc_table.read_from_stream_bin(iif);
    iif.close();
}

//Read cumulative photon emission rate to disk (Warning,
//qs_ctrl is not changed in current implementation)
template<typename _REAL, class _RNDWRAP>
void
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
read_cumulative_phot_em_table(std::string filename)
{
    std::ifstream iif;
    iif.open(filename, std::ios::binary);
    cum_distrib_table.read_from_stream_bin(iif);
    iif.close();
}


//Export innards
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::quantum_synchrotron_innards<_REAL, _RNDWRAP>
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::export_innards()
{
    quantum_synchrotron_innards<_REAL, _RNDWRAP> innards;

    //Filling innards
    innards.qs_ctrl                         = qs_ctrl;
    innards.lambda                          = lambda;
    innards.rng_ptr                         = &rng;
    innards.KKfunc_table_coords_how_many    = KKfunc_table.ref_coords().size();
    innards.KKfunc_table_coords_ptr         = KKfunc_table.ref_coords().data();
    innards.KKfunc_table_data_how_many      = KKfunc_table.ref_data().size();
    innards.KKfunc_table_data_ptr           = KKfunc_table.ref_data().data();
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
picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::quantum_synchrotron_engine
(picsar::multi_physics::quantum_synchrotron_innards<_REAL, _RNDWRAP> innards):
    qs_ctrl{innards.qs_ctrl},
    lambda{innards.lambda},
    KKfunc_table{innards.KKfunc_table_coords_how_many,
                 innards.KKfunc_table_coords_ptr,
                 innards.KKfunc_table_data_ptr},
    cum_distrib_table{innards.cum_distrib_table_coords_1_how_many,
                      innards.cum_distrib_table_coords_1_ptr,
                      innards.cum_distrib_table_coords_2_how_many,
                      innards.cum_distrib_table_coords_2_ptr,
                      innards.cum_distrib_table_data_ptr},
    rng{*innards.rng_ptr}
{}

//___________________PRIVATE FUNCTIONS_____________________________________

//Function to compute y (Warning: it doesn't chek if chi_ele != 0 or
//if chi_phot < chi_ele)
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_y(_REAL chi_phot, _REAL chi_ele) const
{
    return two*chi_phot/(three*chi_ele*(chi_ele-chi_phot));
}

//Function to compute the inner integral of the QS photon emission rate
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_inner_integral(_REAL y) const
{
    auto func = [this](_REAL s){
        if( s > static_cast<_REAL>(__quantum_synchrotron_special_func_big_arg))
            return zero;
        return k_v(five/three, s);
    };

    return quad_a_b<_REAL>(func, y, two*y) + quad_a_inf<_REAL>(func, two*y);
}

template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_KK_integrand(_REAL chi_part, _REAL chi_phot) const
{
    if (chi_part == zero || chi_phot >=  chi_part || chi_phot == zero)
        return zero;

    _REAL y = compute_y(chi_phot, chi_part);

    _REAL inner = compute_inner_integral(y);

    _REAL part_2 = (three*chi_phot*y/two)*k_v(two/three, y);

    _REAL coeff = (one/chi_part)/static_cast<_REAL>(pi*sqrt(three));

    return (inner + part_2)*coeff;
}

template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::quantum_synchrotron_engine<_REAL, _RNDWRAP>::
compute_KK_function(_REAL chi_part) const
{

    auto func = [chi_part, this](_REAL chi_phot){
            return compute_KK_integrand(chi_part, chi_phot);
    };

    //The integral is split into two parts to  handle singularities better
    return quad_a_b<_REAL>(func, zero, chi_part/two) + quad_a_b<_REAL>(func, chi_part/two, chi_part);
}



#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
