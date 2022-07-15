#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CORE

//This .hpp file contains the implementation of the core
//functions of the Quantum Synchrotron photon emission engine.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)


//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Uses picsar arrays
#include "picsar_qed/containers/picsar_array.hpp"
//Uses vector functions
#include "picsar_qed/math/vec_functions.hpp"
//Uses chi functions
#include "picsar_qed/physics/chi_functions.hpp"
//Uses physical constants
#include "picsar_qed/physics/phys_constants.h"
//Uses math constants
#include "picsar_qed/math/math_constants.h"
//Uses unit conversion"
#include "picsar_qed/physics/unit_conversion.hpp"
//Uses sqrt
#include "picsar_qed/math/cmath_overloads.hpp"
//Uses gamma functions
#include "picsar_qed/physics/gamma_functions.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    /**
    * Computes the optical depth of a new electron or positron (simply a number
    * extracted from an exponential distribution).
    *
    * @tparam RealType the floating point type to be used
    * @param[in] unf_zero_one_minus_epsi a random number uniformly distributed in [0,1)
    *
    * @return an exponentially distributed random number
    */
    template<typename RealType>
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    RealType get_optical_depth(const RealType unf_zero_one_minus_epsi)
    {
        using namespace math;
        return -m_log(one<RealType> - unf_zero_one_minus_epsi);
    }

    /**
    * Computes dN/dt for Quantum Synchrotron photon emission. Needs a
    * lookup table to provide G(chi_particle). See validation script
    * for more details.
    *
    * @tparam RealType the floating point type to be used
    * @tparam TableType the type of the lookup table to be used. Must have an "interp" method.
    * @tparam UnitSystem unit system to be used (default is SI)
    *
    * @param[in] t_energy_part particle energy
    * @param[in] chi_part photon chi parameter
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
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    RealType get_dN_dt(const RealType t_energy_part,
        const RealType chi_part,
        const TableType& ref_dndt_table,
        const RealType ref_quantity = math::one<RealType>,
        bool* const is_out_of_table = nullptr)
    {
        const auto energy_part = t_energy_part*conv<
            quantity::energy, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        if(energy_part == math::zero<RealType> || chi_part ==  math::zero<RealType>){
                return  math::zero<RealType>;
        }

        if(chi_part ==  math::zero<RealType>){
                return  math::zero<RealType>;
        }

        const auto GG = ref_dndt_table.interp(chi_part, is_out_of_table);

        constexpr const auto phot_emission_rate_coeff = static_cast<RealType>(
            math::two_thirds<double>*fine_structure<> *
            heaviside_lorentz_electron_rest_energy<double> *
            heaviside_lorentz_electron_rest_energy<double>);

        const auto dndt = phot_emission_rate_coeff * GG/energy_part;

        return dndt*conv<quantity::rate, unit_system::heaviside_lorentz,
            UnitSystem, RealType>::fact(math::one<RealType>,ref_quantity);
    }

    /**
    * Evolves the optical depth of a particle for Quantum Synchrotron pair production.
    * Needs a lookup table to provide G(chi_particle). See validation script
    * for more details.
    *
    * @tparam RealType the floating point type to be used
    * @tparam TableType the type of the lookup table to be used. Must have an "interp" method.
    * @tparam UnitSystem unit system to be used (default is SI)
    *
    * @param[in] t_energy_part particle energy
    * @param[in] chi_part particle chi parameter
    * @param[in] t_dt timestep
    * @param[in,out] the optical depth
    * @param[in] ref_dndt_table a reference to the lookup table
    * @param[in] ref_quantity omega or lambda in SI units if norm_omega or norm_lambda unit systems are used
    *
    * @return true if chi_particle was in the lookup table, false otherwise.
    */
    template<
        typename RealType,
        typename TableType,
        unit_system UnitSystem = unit_system::SI>
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    bool evolve_optical_depth(
        const RealType t_energy_part,
        const RealType chi_part,
        const RealType t_dt, RealType& optical_depth,
        const TableType& ref_dndt_table,
        const RealType ref_quantity = math::one<RealType>)
    {
        const auto energy_part = t_energy_part*conv<
            quantity::energy, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        const auto dt = t_dt*conv<
                quantity::time, UnitSystem,
                unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        bool is_out = false;
        const auto dndt = get_dN_dt<
            RealType, TableType, unit_system::heaviside_lorentz>(
                energy_part, chi_part, ref_dndt_table, ref_quantity,
                &is_out);

        optical_depth -= dndt*dt;

        return !is_out;
    }

    /**
    * Computes the properties of a photon emitted via Quantum Synchrotron photon emission.
    * Updates also the momentum of the emitting particle.
    * Needs a lookup table storing a cumulative probability distribution to
    * calculate the chi parameter of the emitted photon. This lookup table
    * has to provide an "interp" method, accepting the chi of the particle and
    * a uniformly distributed random number in [0,1) as parameters. See validation script
    * for more details.
    *
    * @tparam RealType the floating point type to be used
    * @tparam TableType the type of the lookup table to be used. Must have an "interp" method.
    * @tparam UnitSystem unit system to be used (default is SI)
    *
    * @param[in] chi_particle particle chi parameter
    * @param[in, out] t_v_momentum_particle 3-momentum of the particle
    * @param[in] unf_zero_one_minus_epsi a random number uniformly distributed in [0,1)
    * @param[out] phot_momentum momentum of the generated photon
    * @param[in] ref_quantity omega or lambda in SI units if norm_omega or norm_lambda unit systems are used
    *
    * @return true if chi_particle was in the lookup table, false otherwise.
    */
    template<
        typename RealType, typename TableType,
        unit_system UnitSystem = unit_system::SI
        >
    PXRMP_GPU_QUALIFIER
    PXRMP_FORCE_INLINE
    bool generate_photon_update_momentum(
        const RealType chi_particle,
        math::vec3<RealType>& t_v_momentum_particle,
        const RealType unf_zero_one_minus_epsi,
        const TableType& ref_phot_prod_table,
        math::vec3<RealType>& phot_momentum,
        const RealType ref_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace math;
        using namespace containers;

        const auto mom_u2hl = conv<
            quantity::momentum, UnitSystem,
            unit_system::heaviside_lorentz, RealType>::fact(ref_quantity);

        const auto v_mom_particle = t_v_momentum_particle *mom_u2hl;

        const auto mom_particle = norm(v_mom_particle);
        const auto v_dir_particle = v_mom_particle/mom_particle;
        const auto gamma_particle =
            compute_gamma_ele_pos<RealType, unit_system::heaviside_lorentz>(v_mom_particle);

        bool is_out = false;
        const RealType chi_photon = ref_phot_prod_table.interp(
                chi_particle, unf_zero_one_minus_epsi,
                &is_out);

        if (chi_photon == zero<RealType>){
            phot_momentum = math::vec3<RealType>{
                zero<RealType>, zero<RealType>, zero<RealType>};
            return false;
        }

        const auto gamma_photon = (gamma_particle - one<RealType>)*
            chi_photon/chi_particle;

        const auto mom_hl2u = conv<
            quantity::momentum, unit_system::heaviside_lorentz,
            UnitSystem, RealType>::fact(one<RealType>, ref_quantity);

        phot_momentum = heaviside_lorentz_electron_rest_energy<RealType>*
            v_dir_particle*gamma_photon*mom_hl2u;
        t_v_momentum_particle = t_v_momentum_particle - phot_momentum;

        return !is_out;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CORE
