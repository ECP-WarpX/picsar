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

    //______________________GPU
    //This function computes the properties of the electron-positron pairs
    //generated in a BW process.
    //Conceived for GPU usage.
    template<
        typename RealType, typename TableType,
        unit_system UnitSystem = unit_system::SI
        >
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    void generate_breit_wheeler_pairs(
        const RealType chi_photon,
        const math::vec3<RealType>& t_v_momentum_photon,
        const RealType unf_zero_one_minus_epsi,
        const TableType& ref_pair_prod_table,
        math::vec3<RealType>& ele_momentum,
        math::vec3<RealType>& pos_momentum,
        const RealType ref_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace math;
        using namespace containers;

        const auto mon_u2hl = conv<
            quantity::momentum, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        const auto v_mom_photon = t_v_momentum_photon *mon_u2hl;

        const auto mom_photon = norm(v_mom_photon);
        const auto v_dir_photon = v_mom_photon/mom_photon;
        const auto gamma_photon = mom_photon/
            heaviside_lorentz_electron_rest_energy<RealType>;

        const auto chi_ele = ref_pair_prod_table.interp(
                chi_photon, unf_zero_one_minus_epsi);
        const auto chi_pos = chi_photon - chi_ele;

        const auto coeff =  (gamma_photon - two<RealType>)/chi_photon;
        const auto gamma_ele = one<RealType>+chi_ele*coeff;
        const auto gamma_pos = one<RealType>+chi_pos*coeff;

        const auto mom_hl2u = conv<
            quantity::momentum, unit_system::heaviside_lorentz,
            UnitSystem, RealType>::fact(one<RealType>, ref_quantity);

        ele_momentum = heaviside_lorentz_electron_rest_energy<RealType>*
            v_dir_photon*sqrt(gamma_ele*gamma_ele - one<RealType>)*mom_hl2u;
        pos_momentum = heaviside_lorentz_electron_rest_energy<RealType>*
            v_dir_photon*sqrt(gamma_pos*gamma_pos - one<RealType>)*mom_hl2u;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
