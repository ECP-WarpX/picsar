#ifndef PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
#define PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION

#include "../qed_commons.h"
#include "phys_constants.h"
#include "../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace phys{
    enum unit_system
    {
        SI,
        norm_omega,
        norm_lambda
    };

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_mass_to_SI_from(){
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::electron_mass);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(phys::electron_mass);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_mass_from_SI_to(){
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(1.0/phys::electron_mass);
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0/phys::electron_mass);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_charge_to_SI_from(){
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::elementary_charge);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(phys::elementary_charge);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_charge_from_SI_to(){
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(1.0/phys::elementary_charge);
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0/phys::elementary_charge);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_velocity_to_SI_from(){
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::light_speed);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_velocity_from_SI_to(){
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(1.0/phys::light_speed);
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0/phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

        template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_momentum_to_SI_from(){
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::electron_mass*phys::light_speed);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(phys::electron_mass*phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_momentum_from_SI_to(){
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(1.0/(phys::electron_mass*phys::light_speed));
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0/(phys::electron_mass*phys::light_speed));
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_length_to_SI_from(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return lambda*static_cast<RealType>(1.0/2.0/math::pi);
        else if (From == unit_system::norm_lambda)
            return lambda;
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_length_from_SI_to(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(2.0*math::pi)/lambda;
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0)/lambda;
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_time_to_SI_from(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return lambda*static_cast<RealType>(1.0/2.0/phys::light_speed/math::pi);
        else if (From == unit_system::norm_lambda)
            return lambda*static_cast<RealType>(1.0/phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_time_from_SI_to(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(2.0*math::pi*phys::light_speed)/lambda;
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(phys::light_speed)/lambda;
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_rate_to_SI_from(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(2.0*phys::light_speed*math::pi)/lambda;
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(phys::light_speed)/lambda;
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_rate_from_SI_to(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return lambda*static_cast<RealType>(1.0/2.0/math::pi/phys::light_speed);
        else if (To == unit_system::norm_lambda)
            return lambda*static_cast<RealType>(1.0/phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_E_to_SI_from(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(
                (phys::electron_mass*phys::light_speed*2.0*math::pi*phys::light_speed)/
                (phys::elementary_charge))/lambda;
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(
                phys::electron_mass*phys::light_speed*phys::light_speed/
                (phys::elementary_charge))/lambda;
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_E_from_SI_to(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return lambda*static_cast<RealType>(
                (phys::elementary_charge)/
                (phys::electron_mass*phys::light_speed*2.0*math::pi*phys::light_speed));
        else if (To == unit_system::norm_lambda)
            return lambda*static_cast<RealType>(
                (phys::elementary_charge)/
                (phys::electron_mass*phys::light_speed*phys::light_speed));
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_B_to_SI_from(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(
                (phys::electron_mass*2.0*math::pi*phys::light_speed)/
                (lambda*phys::elementary_charge));
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(
                phys::electron_mass*phys::light_speed/
                (lambda*phys::elementary_charge));
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_B_from_SI_to(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(
                (lambda*phys::elementary_charge)/
                (phys::electron_mass*2.0*math::pi*phys::light_speed));
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(
                (lambda*phys::elementary_charge)/
                (phys::electron_mass*phys::light_speed));
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_energy_to_SI_from()
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(
                phys::electron_mass*phys::light_speed*phys::light_speed);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(
                phys::electron_mass*phys::light_speed*phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_energy_from_SI_to(
        const RealType lambda = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(
                1.0/phys::electron_mass*phys::light_speed*phys::light_speed);
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(
                1.0/phys::electron_mass*phys::light_speed*phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

}
}
}
#endif //PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION