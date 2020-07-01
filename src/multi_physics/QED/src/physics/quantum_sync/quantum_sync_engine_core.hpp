#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CORE
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CORE

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

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

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
    RealType get_dN_dt(const RealType t_energy_part,
        const RealType chi_part,
        const TableType& ref_dndt_table,
        const RealType ref_quantity = math::one<RealType>)
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

        const auto GG = ref_dndt_table.interp(chi_part);

        constexpr const auto phot_emission_rate_coeff = static_cast<RealType>(
            math::two_thirds<double>*fine_structure<> *
            heaviside_lorentz_electron_rest_energy<double> *
            heaviside_lorentz_electron_rest_energy<double>);

        const auto dndt = phot_emission_rate_coeff * GG/energy_part;

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

        const auto dndt = get_dN_dt<
            RealType, TableType, unit_system::heaviside_lorentz>(
                energy_part, chi_part, ref_dndt_table, ref_quantity);

        optical_depth -= dndt*dt;

        return (optical_depth <= math::zero<RealType>);
    }



}
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_CORE
