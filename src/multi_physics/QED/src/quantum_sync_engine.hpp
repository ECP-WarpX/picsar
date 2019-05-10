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
          size_t chi_frac_tem_how_many =
            __quantum_synchrotron_chi_frac_tem_how_many;
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
        lambda*one/(energy_part*chi_part);

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
    picsar_vector<_REAL> KK_vals{};
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
        _lambda*one/(energy_part*chi_part);

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
    auto func = [this, chi_part](_REAL chi_phot){
        return compute_KK_integrand(chi_phot, chi_part);
    };
   _REAL num = quad_a_b<_REAL>(func, zero, chi_phot);
   return num/compute_KK_function(chi_part) ;
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

    picsar_vector<_REAL> frac_coords = generate_lin_spaced_vec(zero, one,
    qs_ctrl.chi_frac_tem_how_many);

    picsar_vector<_REAL> pair_vals{};

    msg("Computing table for photon emission...\n", stream);

    for(auto chi_part: chi_coords){
        pair_vals.push_back(zero);
        msg("chi_part: " + std::to_string(chi_part) + " \n", stream);
        for(size_t i = 1; i < frac_coords.size() - 1; i++){
            _REAL temp = compute_cumulative_phot_em(
                chi_part*frac_coords[i], chi_part);
            pair_vals.push_back(temp);

        }
        pair_vals.push_back(one); //The function is symmetric
    }
    msg("...done!\n", stream);

    cum_distrib_table = lookup_2d<_REAL>{
        picsar_array<picsar_vector<_REAL>,2>{chi_coords, frac_coords}, pair_vals};


    //Initialize the auxiliary table
    aux_table = lookup_1d<_REAL>{cum_distrib_table.get_coords()[1],
    picsar_vector<_REAL>(qs_ctrl.chi_frac_tem_how_many)};
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
    picsar_vector<std::pair<vec3<_REAL>, _REAL>> photons(sampling);

    _REAL chi_part = chi_lepton(px, py, pz, ex, ey, ez, bx, by, bz, lambda);

    auto frac = aux_table.ref_coords();

    //if chi < chi_min: chi = chi_min for cumulative distribution
    if(chi_part < qs_ctrl.chi_part_tem_min){
        for(size_t i = 0; i < frac.size(); i++){
            aux_table.ref_data()[i] =
            cum_distrib_table.data_at_coords(0, i);
        }
    }
    //if chi > chi_max: chi = chi_max for cumulative distribution
    else if(chi_part < qs_ctrl.chi_part_tem_max){
        for(size_t i = 0; i < frac.size(); i++){
            aux_table.ref_data()[i] =
            cum_distrib_table.data_at_coords(qs_ctrl.chi_part_tem_how_many, i);
        }
    }
    //Interpolate 1D cumulative distribution
    else{
        for(size_t i = 0; i < frac.size(); i++){
            aux_table.ref_data()[i] =
            cum_distrib_table.interp_linear_first_equispaced(chi_part, i);
        }
    }

    _REAL new_weight = weight/sampling;

    _REAL me_c = static_cast<_REAL>(__emass*__c);

    vec3<_REAL> p_part{px, py, pz};
    _REAL norm_part = norm(p_part);
    vec3<_REAL> n_part = p_part/norm_part;
    _REAL gamma_part = norm_part/me_c;

    _REAL one_over_sampling = one/sampling;

    for(size_t s = 0; s < sampling; s++){
        _REAL prob = rng.unf(zero, one);

        _REAL chi_phot = aux_table.interp_linear_equispaced(prob);

        _REAL gamma_phot = chi_phot/chi_part*(gamma_part-one);

        vec3<_REAL>  p_phot = gamma_phot * n_part * me_c;

        photons[s] = std::make_pair(p_phot, new_weight);

        //Update particle momentum
        px -= p_phot[0]*one_over_sampling;
        py -= p_phot[1]*one_over_sampling;
        pz -= p_phot[2]*one_over_sampling;
    }

    return photons;
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

    aux_table = lookup_1d<_REAL>{cum_distrib_table.get_coords()[1],
    picsar_vector<_REAL>(qs_ctrl.chi_frac_tem_how_many)};
}


//___________________PRIVATE FUNCTIONS_____________________________________

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
    if (chi_part == zero || chi_phot == zero)
        return zero;

    _REAL y = compute_y(chi_phot, chi_part);

    _REAL inner = compute_inner_integral(y);

    _REAL part_2 = (two + three*chi_phot*y)*k_v(two/three, y);

    _REAL coeff = one/static_cast<_REAL>(pi*sqrt(three));

    return -(inner - part_2)*coeff;
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

    return quad_a_b<_REAL>(func, zero, chi_part);
}



#endif //__PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
