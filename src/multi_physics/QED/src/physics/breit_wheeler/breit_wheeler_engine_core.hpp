#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CORE

//This .hpp file contais the implementation of the
//core functions of the nonlinear Breit Wheeler engine.
//If desired, these functions can be used directly
//without the higher-level interface.

#include <cmath>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar arrays
#include "../../containers/picsar_array.hpp"

//Uses vector functions
#include "../../math/vec_functions.hpp"

//Uses chi functions
#include "../chi_functions.hpp"

//Uses physical constants
#include "../phys_constants.h"

//Uses unit conversion"
#include "../unit_conversion.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    enum check_input
    {
        yes,
        no
    };

    template<typename RealType>
    struct control_params{
        //Minimum chi_phot to consider in all the functions
        //If chi is less than __bw_min_chi_phot
        //the function evolve_opt_depth_and_determine_event
        //provided by the BW engine will return immediately
        //without touching the optical depth and communicating that no event
        //has occurred.
        RealType chi_phot_min = static_cast<RealType>(0.001);
    };

    template<typename RealType>
    struct dndt_lookup_table_params{
        // dN/dt table:
        //__breit_wheeler_min_tdndt_chi_phot  is the inferior
        //limit of the total pair production rate lookup table.
        //If  __breit_wheeler_min_chi_phot < chi <__breit_wheeler_min_tdndt_chi_phot
        //BW process is taken into account, but the Erber approximation
        //rather than the lookup table is used.
        RealType chi_phot_tdndt_min = static_cast<RealType>(0.1); //Min chi_phot
        RealType chi_phot_tdndt_max = static_cast<RealType>(500.0); //Max chi_phot
        size_t chi_phot_tdndt_how_many = 50; //How many points
    };


    template<typename RealType, typename VectorType>
    class dndt_lookup_table{
        public:
            template <unit_system UnitSystem, check_input CheckInput>
            friend RealType interp_dN_dt(
                const RealType, const RealType,
                const dndt_lookup_table<RealType, VectorType>&,
                const control_params<RealType>&);

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(RealType where) const
            {
                const auto xmin = m_ctrl.chi_phot_tdndt_min;
                const auto xmax = m_ctrl.chi_phot_tdndt_max;
                const auto isize = m_ctrl.chi_phot_tdndt_how_many;
                const auto xsize = xmax - xmin;

                const size_t idx_left = static_cast<size_t>(
                    floor((isize-1)*(where-xmin)/xsize));
                const size_t idx_right = idx_left + 1;

                const auto xleft = (idx_left*xsize/isize) + xmin;
                const auto xright = (idx_right*xsize/isize) + xmin;
                const auto yleft = m_values[idx_left];
                const auto yright = m_values[idx_right];

                return yleft + ((yright-yleft)/(xright-xleft))*(where-xleft);
            }

        private:
            dndt_lookup_table_params<RealType> m_ctrl;
            VectorType m_values;
    };


    template <typename T, typename B>
    T funct (T a, B c){return a + c;}

    //______________________GPU
    //get a single optical depth
    //same as above but conceived for GPU usage
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_optical_depth(const RealType unf_zero_one_minus_epsi)
    {
        return -log(static_cast<RealType>(1.0) - unf_zero_one_minus_epsi);
    }

    //______________________GPU
    //Interp the dN_dt from table (with GPU use directly this one!)
    template<
        typename RealType,
        typename VectorType,
        check_input CheckInput = check_input::yes,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType interp_dN_dt(
        const RealType t_energy_phot, const RealType chi_phot,
        const dndt_lookup_table<RealType, VectorType>& ref_dndt_table,
        const control_params<RealType>& ref_ctrl,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        constexpr const auto zero = static_cast<RealType>(0.0);
        constexpr const auto one = static_cast<RealType>(1.0);
        constexpr const auto two = static_cast<RealType>(2.0);
        constexpr const auto three = static_cast<RealType>(3.0);

        const auto energy_phot = t_energy_phot*
            fact_energy_to_SI_from<UnitSystem, RealType>(ref_quantity);

        PXRMP_INTERNAL_CONSTEXPR_IF (CheckInput == check_input::yes){
            if(energy_phot == zero || chi_phot == zero)
                return zero;
        }

        constexpr const auto pair_prod_rate_coeff = static_cast<RealType>(
            phys::fine_structure * phys::electron_mass *
            phys::light_speed * phys::light_speed / reduced_plank);

        const auto coeff = pair_prod_rate_coeff/( chi_phot * energy_phot);

        const auto TT = zero;
        //Use approximate analytical expression if chi < chi_phot_tdndt_min
        //or chi > chi_phot_tdndt_min

        //Coefficients for che asymptotic behaviour of the "TT function"
        //using Erber approximation
        //Erber T. (1966), Reviews of Modern Physics,38, 626
        constexpr const auto erber_a = static_cast<RealType>(0.16);
        constexpr const auto erber_b = static_cast<RealType>(4./3.);

        if(chi_phot <= ref_dndt_table.ctrl.chi_phot_tdndt_min){
            TT = (static_cast<RealType>(math::pi)*erber_a/erber_b/two)*
                chi_phot*chi_phot*exp(-two*erber_b/chi_phot);
        }
        else if(chi_phot >= ref_dndt_table.ctrl.chi_phot_tdndt_max){
            constexpr const auto ccoeff =
                static_cast<RealType>(tgamma(one/three)/two);
            TT = erber_a*chi_phot*ccoeff*pow(chi_phot*two/erber_b, two/three);
        }
        //otherwise use lookup tables
        else{
            TT =  exp(ref_dndt_table.interp(log(chi_phot)));
        }

        const auto dndt = coeff * TT;

        return dndt*fact_rate_from_SI_to<UnitSystem, RealType>(ref_quantity);
    }

    //______________________GPU
    //This function evolves the optical depth for a particle and
    //checks if it goes to zero.
    //Conceived for GPU usage.
    template<
        typename RealType,
        typename VectorType,
        check_input CheckInput = check_input::yes,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    void evolve_opt_depth_and_determine_event(
        const RealType chi_photon, const RealType t_energy_photon,
        const RealType t_dt, RealType& opt_depth,
        bool& has_event_happend, RealType& event_dt,
        const dndt_lookup_table<RealType, VectorType>& ref_dndt_table,
        const control_params<RealType>& ref_ctrl,
        const RealType ref_quantity = static_cast<RealType>(1.0))
    {
        constexpr const auto zero = static_cast<RealType>(0.0);
        constexpr const auto two = static_cast<RealType>(2.0);

        const auto energy_photon = t_energy_photon*
            fact_energy_to_SI_from<UnitSystem, RealType>(ref_quantity);

        const auto dt = t_dt*
            fact_time_to_SI_from<UnitSystem, RealType>(ref_quantity);

        has_event_happend = false;
        event_dt = zero;

        //Do not evolve opt_depth if the photon energy is not high enough to generate a pair
        PXRMP_INTERNAL_CONSTEXPR_IF (CheckInput == check_input::yes){
            constexpr const auto energy_threshold =
                two*phys::electron_mass*phys::light_speed*phys::light_speed;
            if(t_energy_photon < energy_threshold)
                return;
        }

        //Do not evolve opt_depth if the chi parameter is less then threshold
        if(chi_photon <= ref_ctrl.chi_phot_min)
            return;

        //**Compute dndt (internally SI units are used)
        const auto dndt = interp_dN_dt<RealType, VectorType, CheckInput, unit_system::SI>(
                energy_photon, chi_photon, ref_dndt_table, ref_ctrl, ref_quantity);

        opt_depth -= dndt*dt;

        if(opt_depth < zero){
            //Calculates the time at which pair prod takes place
            event_dt = (opt_depth/dndt + dt)*
                fact_time_from_SI_to<UnitSystem, RealType>(ref_quantity);
            has_event_happend = true;
        }
    }

    //______________________GPU
    //This function computes the properties of the electron-positron pairs
    //generated in a BW process.
    //Conceived for GPU usage.
    template<
        typename RealType,
        typename VectorType,
        check_input CheckInput = check_input::yes,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    void generate_breit_wheeler_pairs(
        const RealType chi_photon, const RealType t_momentum_photon,
        const size_t sampling, RealType* ele_momentum, RealType* pos_momentum,
        const RealType* unf_zero_one_minus_epsi,
        const dndt_lookup_table<RealType, VectorType>& ref_dndt_table,
        const control_params<RealType>& ref_ctrl,
        const RealType ref_quantity = static_cast<RealType>(1.0))
{
    constexpr const auto me_c =
        static_cast<RealType>(phys::electron_mass*phys::light_speed);

    const auto momentum_photon = t_momentum_photon*
        fact_momentum_to_SI_from<UnitSystem,RealType>();

    const auto gamma_phot = momentum_photon/me_c;

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
    }
}


}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
