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

//Uses math constants
#include "../../math/math_constants.h"

//Uses unit conversion"
#include "../unit_conversion.hpp"

#include "breit_wheeler_engine_tables.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    //______________________GPU
    //get a single optical depth
    //same as above but conceived for GPU usage
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_optical_depth(const RealType unf_zero_one_minus_epsi)
    {
        using namespace math;
        return -log(one<RealType> - unf_zero_one_minus_epsi);
    }


    //______________________GPU
    //Interp the dN_dt from table (with GPU use directly this one!)
    template<
        typename RealType,
        typename TableType,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_dN_dt(
        const RealType t_energy_phot, const RealType chi_phot,
        const TableType& ref_dndt_table,
        const RealType ref_quantity = math::one<RealType>)
    {
        const auto energy_phot = t_energy_phot*conv<
            quantity::energy, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        if(energy_phot == math::zero<RealType> || chi_phot ==  math::zero<RealType>){
                return  math::zero<RealType>;
        }

        const auto TT = ref_dndt_table.interp(chi_phot);

        constexpr const auto pair_prod_rate_coeff = static_cast<RealType>(
            fine_structure<> * heaviside_lorentz_electron_rest_energy<RealType>
            * heaviside_lorentz_electron_rest_energy<RealType>);

        const auto dndt = pair_prod_rate_coeff * TT *(chi_phot/energy_phot);

        return dndt*conv<quantity::rate, unit_system::heaviside_lorentz,
            UnitSystem, RealType>::fact(math::one<RealType>,ref_quantity);
    }

    //______________________GPU
    //Evolve optical depth (with GPU use directly this one!)
    template<
        typename RealType,
        typename TableType,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    bool evolve_optical_depth(
        const RealType t_energy_phot, const RealType chi_phot,
        const RealType t_dt, RealType& optical_depth,
        const TableType& ref_dndt_table,
        const RealType ref_quantity = math::one<RealType>)
    {
        const auto energy_phot = t_energy_phot*conv<
            quantity::energy, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        const auto dt = t_dt*conv<
                quantity::time, UnitSystem,
                unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        const auto dndt = get_dN_dt<
            RealType, TableType, unit_system::heaviside_lorentz>(
                energy_phot, chi_phot, ref_dndt_table, ref_quantity);

        optical_depth -= dndt*dt;

        return (optical_depth <= math::zero<RealType>);
    }

/*

    template<typename RealType>
    struct control_params{
        //Minimum chi_phot to consider in all the functions
        //If chi is less than __bw_min_chi_phot
        //the function evolve_opt_depth_and_determine_event
        //provided by the BW engine will return immediately
        //without touching the optical depth and communicating that no event
        //has occurred.
        RealType chi_phot_min = static_cast<RealType>(0.001);

        //Sets the limits for the semi-infinite integrals in the library
        RealType very_large_number = 1.0e20;
    };



    template<typename RealType>
    struct pair_prod_lookup_table_params{
        // Pair production table:
        RealType chi_phot_tpair_min = static_cast<RealType>(0.01); //Min chi_phot
        RealType chi_phot_tpair_max = static_cast<RealType>(500.0); //Max chi_phot
        size_t chi_phot_tpair_how_many = 50; //How many points
        size_t chi_frac_tpair_how_many = 50;
    };


    template<typename RealType, typename VectorType>
    class pair_prod_lookup_table{
        public:
            template<unit_system UnitSystem>
            friend void generate_breit_wheeler_pairs(
                const RealType, const RealType,
                const size_t, RealType*, RealType*,
                const RealType*, const pair_prod_lookup_table<RealType, VectorType>&,
                const control_params<RealType>&, const RealType);

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp_linear_first(RealType where_x, size_t coord_y) const
            {
                const auto sx = m_ctrl.chi_phot_tpair_how_many;
                const auto sy = m_ctrl.chi_frac_tpair_how_many;

                const auto xmin = m_ctrl.chi_phot_tpair_min;
                const auto xmax = m_ctrl.chi_phot_tpair_max;
                const auto xsize = xmax - xmin;

                if(where_x <= xmin)
                    return m_values[idx(0,coord_y)];

                if(where_x >= xmax)
                    return m_values[idx(sx-1,coord_y)];

                const auto idx_x_left = static_cast<size_t>(
                    floor((sx-1)*(where_x-xmin)/xsize));
                const auto idx_x_right = idx_x_left + 1;

                const auto xleft = (idx_x_left*xsize/sx) + xmin;
                const auto xright = (idx_x_right*xsize/sx) + xmin;

                const auto zlc = m_values[idx(idx_x_left,coord_y)];
                const auto zrc = m_values[idx(idx_x_right,coord_y)];

                const auto wlc = (xright - where_x);
                const auto wrc = (where_x - xleft);

                const auto w_norm = (xright-xleft);

                return (wlc*zlc + wrc*zrc)/w_norm;

            }


        private:
            pair_prod_lookup_table_params<RealType> m_ctrl;
            VectorType m_values;

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            int idx(int i, int j) const
            {
                return i*m_ctrl.chi_phot_tpair_how_many + j;
            };
    };




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
        typename RealType, typename VectorType,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    void generate_breit_wheeler_pairs(
        const RealType chi_photon, const RealType t_momentum_photon,
        const size_t sampling,
        PXRMR_INTERNAL_RESTRICT RealType* ele_momentum,
        PXRMR_INTERNAL_RESTRICT RealType* pos_momentum,
        PXRMR_INTERNAL_RESTRICT const RealType* unf_zero_one_minus_epsi,
        const pair_prod_lookup_table<RealType, VectorType>& ref_pair_prod_table,
        const control_params<RealType>& ref_ctrl,
        const RealType ref_quantity = static_cast<RealType>(1.0))
{
    constexpr const auto zero = static_cast<RealType>(0.0);
    constexpr const auto one = static_cast<RealType>(1.0);
    constexpr const auto two = static_cast<RealType>(2.0);

    constexpr const auto me_c =
        static_cast<RealType>(phys::electron_mass*phys::light_speed);

    const auto momentum_photon = t_momentum_photon*
        fact_momentum_to_SI_from<UnitSystem,RealType>();

    const auto gamma_phot = momentum_photon/me_c;

    const size_t how_many_frac =
        ref_pair_prod_table.m_ctrl.chi_frac_tpair_how_many;

    auto tab_chi_phot = chi_photon;
    if(chi_photon < ref_pair_prod_table.m_ctrl.chi_phot_tpair_min)
        tab_chi_phot = ref_pair_prod_table.m_ctrl.chi_phot_tpair_min;
    else if(chi_photon > ref_pair_prod_table.m_ctrl.chi_phot_tpair_max)
        tab_chi_phot = ref_pair_prod_table.m_ctrl.chi_phot_tpair_max;

    for(size_t s = 0; s < sampling; ++s){
        auto prob = unf_zero_one_minus_epsi[s];
        bool invert = false;
        if(prob >= one/two){
            prob -= one/two;
            invert = true;
        }

        size_t first = 0;
        auto val = zero;
        size_t count = how_many_frac;
        while(count > 0){
            auto it = first;
            size_t step = count/2;
            it += step;
            val =
                exp(ref_pair_prod_table.interp_linear_first
                    (log(tab_chi_phot), it));
            if(!(prob < val)){
                first = ++it;
                count -= step+1;
            }
            else{
                count = step;
            }
        }

        const auto upper = first;
        const auto lower = upper-1;

        const auto upper_frac =
            upper*0.5/(ref_pair_prod_table.m_ctrl.chi_frac_tpair_how_many-1);

        const auto lower_frac =
            upper*0.5/(ref_pair_prod_table.m_ctrl.chi_frac_tpair_how_many-1);

        const auto upper_prob =
            exp(ref_pair_prod_table.interp_linear_first(
                log(tab_chi_phot), upper));
        const auto lower_prob =
            exp(ref_pair_prod_table.interp_linear_first(
                log(tab_chi_phot), lower));
        const auto chi_ele_frac = lower_frac +
            (prob-lower_prob)*(upper_frac-lower_frac)/(upper_prob-lower_prob);

        const auto chi_ele = chi_ele_frac*chi_photon;
        const auto chi_pos = chi_photon - chi_ele;

        if(invert){
            const auto tt = chi_ele;
            chi_ele = chi_pos;
            chi_pos = tt;
        }

        const auto coeff =  (gamma_phot - two)/chi_photon;

        ele_momentum[s] = sqrt((one+chi_ele*coeff)*(one+chi_ele*coeff)-one)*me_c*
            fact_momentum_from_SI_to<UnitSystem,RealType>();
        pos_momentum[s] = sqrt((one+chi_pos*coeff)*(one+chi_pos*coeff)-one)*me_c*
            fact_momentum_from_SI_to<UnitSystem,RealType>();;
    }
}
*/

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
