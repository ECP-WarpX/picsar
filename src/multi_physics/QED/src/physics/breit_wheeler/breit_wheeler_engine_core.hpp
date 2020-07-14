#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_CORE

//This .hpp file contais the implementation of the core
//function of the Breit-Wheeler pair production engine.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) A.I.Nikishov. & V.I. Ritus Sov. Phys. JETP 19, 2 (1964)
// 2) T.Erber Rev. Mod. Phys. 38, 626 (1966)
// 3) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 4) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses GPU-friendly arrays
#include "../../containers/picsar_array.hpp"
//Uses vector functions
#include "../../math/vec_functions.hpp"
//Uses chi functions
#include "../chi_functions.hpp"
//Uses physical constants
#include "../phys_constants.h"
//Uses mathematical constants
#include "../../math/math_constants.h"
//Uses unit conversion
#include "../unit_conversion.hpp"
//Uses sqrt and log
#include "../../math/cmath_overloads.hpp"

#include <cmath>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    /**
    * Computes the optical depth of a new photon (simply a number
    * extracted from an exponential distribution).
    *
    * @tparam RealType the floating point type to be used
    * @param[in] unf_zero_one_minus_epsi a random number uniformly distributed in [0,1)
    *
    * @return an exponentially distributed random number
    */
    template<typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_optical_depth(const RealType unf_zero_one_minus_epsi)
    {
        using namespace math;
        return -m_log(one<RealType> - unf_zero_one_minus_epsi);
    }


    /**
    * Computes dN/dt for Breit-Wheeler pair production. Needs a
    * lookup table to provide T(chi_photon). See validation script
    * for more details.
    *
    * @tparam RealType the floating point type to be used
    * @tparam TableType the type of the lookup table to be used. Must have an "interp" method.
    * @tparam UnitSystem unit system to be used (default is SI)
    *
    * @param[in] t_energy_phot photon energy
    * @param[in] chi_phot photon chi parameter
    * @param[in] ref_dndt_table a reference to the lookup table
    * @param[in] ref_quantity omega or lambda in SI units if norm_omega or norm_lambda unit systems are used
    * @param[out] is_out_of_table if provided it is set to true in case chi_phot is out of table
    *
    * @return total pair production cross section dN/dt in UnitSystem
    */
    template<
        typename RealType,
        typename TableType,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType get_dN_dt(
        const RealType t_energy_phot, const RealType chi_phot,
        const TableType& ref_dndt_table,
        const RealType ref_quantity = math::one<RealType>,
        bool* const is_out_of_table = nullptr)
    {
        const auto energy_phot = t_energy_phot*conv<
            quantity::energy, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        if(energy_phot == math::zero<RealType> || chi_phot ==  math::zero<RealType>){
                return  math::zero<RealType>;
        }

        const auto TT = ref_dndt_table.interp(chi_phot, is_out_of_table);

        constexpr const auto pair_prod_rate_coeff = static_cast<RealType>(
            fine_structure<> * phys::heaviside_lorentz_electron_rest_energy<RealType>
            * phys::heaviside_lorentz_electron_rest_energy<RealType>);

        const auto dndt = pair_prod_rate_coeff * TT *(chi_phot/energy_phot);

        return dndt*conv<quantity::rate, unit_system::heaviside_lorentz,
            UnitSystem, RealType>::fact(math::one<RealType>,ref_quantity);
    }

    /**
    * Evolves the optical depth of a photon for Breit-Wheeler pair production.
    * Needs a lookup table to provide T(chi_photon). See validation script
    * for more details.
    *
    * @tparam RealType the floating point type to be used
    * @tparam TableType the type of the lookup table to be used. Must have an "interp" method.
    * @tparam UnitSystem unit system to be used (default is SI)
    *
    * @param[in] t_energy_phot photon energy
    * @param[in] chi_phot photon chi parameter
    * @param[in] t_dt timestep
    * @param[in,out] the optical depth
    * @param[in] ref_dndt_table a reference to the lookup table
    * @param[in] ref_quantity omega or lambda in SI units if norm_omega or norm_lambda unit systems are used
    *
    * @return true if chi_phot was in the lookup table, false otherwise.
    */
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

        bool is_out = false;
        const auto dndt = get_dN_dt<
            RealType, TableType, unit_system::heaviside_lorentz>(
                energy_phot, chi_phot, ref_dndt_table, ref_quantity,
                &is_out);

        optical_depth -= dndt*dt;

        return !is_out;
    }

    /**
    * Computes the properties of particles emitted in Breit-Wheeler pair production.
    * Needs a lookup table storing a cumulative probability distribution to
    * calculate the chi parameter of the emitted electron. This lookup table
    * has to provide an "interp" method, accepting the chi of the photon and
    * a uniformly distributed random number in [0,1) as parameters. See validation script
    * for more details.
    *
    * @tparam RealType the floating point type to be used
    * @tparam TableType the type of the lookup table to be used. Must have an "interp" method.
    * @tparam UnitSystem unit system to be used (default is SI)
    *
    * @param[in] chi_phot photon chi parameter
    * @param[in] t_v_momentum_photon 3-momentum of the photon
    * @param[in] unf_zero_one_minus_epsi a random number uniformly distributed in [0,1)
    * @param[out] ele_momentum momentum of the generated electron
    * @param[out] pos_momentum momentum of the generated positron
    * @param[in] ref_quantity omega or lambda in SI units if norm_omega or norm_lambda unit systems are used
    *
    * @return true if chi_phot was in the lookup table, false otherwise.
    */
    template<
        typename RealType, typename TableType,
        unit_system UnitSystem = unit_system::SI
        >
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    bool generate_breit_wheeler_pairs(
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

        bool is_out = false;
        const auto chi_ele = ref_pair_prod_table.interp(
                chi_photon, unf_zero_one_minus_epsi,
                &is_out);
        const auto chi_pos = chi_photon - chi_ele;

        const auto coeff =  (gamma_photon - two<RealType>)/chi_photon;
        const auto gamma_ele = one<RealType>+chi_ele*coeff;
        const auto gamma_pos = one<RealType>+chi_pos*coeff;

        const auto mom_hl2u = conv<
            quantity::momentum, unit_system::heaviside_lorentz,
            UnitSystem, RealType>::fact(one<RealType>, ref_quantity);

        ele_momentum = phys::heaviside_lorentz_electron_rest_energy<RealType>*
            v_dir_photon*m_sqrt(gamma_ele*gamma_ele - one<RealType>)*mom_hl2u;
        pos_momentum = phys::heaviside_lorentz_electron_rest_energy<RealType>*
            v_dir_photon*m_sqrt(gamma_pos*gamma_pos - one<RealType>)*mom_hl2u;

        return !is_out;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE_CORE
