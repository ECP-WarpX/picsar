#ifndef PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
#define PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION

#include "../qed_commons.h"
#include "phys_constants.h"
#include "../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace phys{

    /**
    * The interface of the PICSAR QED library supports 4 unit systems
    * (Internally all the calculations are performed in natural units):
    *
    * > SI units (International System of Units)
    *
    * > Natural units, a unit system frequently used in particle phyics
    *   where c (speed of light) = hbar (reduced Plank constant) = 1.
    *   We also choose to measure energy in GeV (1 GeV ~ 1.60218e-10 J ) and
    *   to use Heaviside-Lorentz units for charge and electromagnetic fields.    *
    *   As a consequence:
    *   - mass is in 1/GeV units
    *   - e
    *
    * > norm_omega, a popular unit system in Particle-In-Cell codes where:
    *   - a reference frequenecy omega_r = 2*pi*c/lambda is chosen
    *     (c is the speed of light, lambda is a wavelength)
    *   - mass is normalized with respect to the electron mass m_e
    *   - electric charge is normalized with respect to the elementary charge e0
    *   - momentum is normalized with respect to m_e*c
    *   - energy is normalized with respect to m_e*c^2
    *   - velocity is normalized with respect to the speed of light
    *   - time is normalized with respect to 1/omega_r
    *   - length is normalized with respect to c/omega_r
    *   - electric field is normalized with respect to m_e*c*omega_r /e0
    *   - magnetic field is normalized with respect to m_e*omega_r /e0
    *
    * > norm_lambda, a popular unit system in Particle-In-Cell codes where:
    *   - a reference wavelength lambda is chosen
    *   - mass is normalized with respect to the electron mass m_e
    *   - electric charge is normalized with respect to the elementary charge e0
    *   - momentum is normalized with respect to m_e*c
    *   - energy is normalized with respect to m_e*c^2
    *   - velocity is normalized with respect to the speed of light
    *   - time is normalized with respect to lambda/c
    *   - length is normalized with respect to lambda
    *   - electric field is normalized with respect to m_e*c^2/(lambda*e0)
    *   - magnetic field is normalized with respect to m_e*c/(lambda*e0)
    *
    * The helper functions provided by this module are used internally by the
    * library to convert the input values of a function to SI and to convert the
    * output from SI to the desired unit system.
    *
    */
    enum unit_system
    {
        SI,
        natural,
        norm_omega,
        norm_lambda
    };

    /**
    * This function returns the conversion factor for mass from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for mass from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for charge from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for charge from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for velocity from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for velocity from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for momentum from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for momentum from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for length from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_length_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::light_speed /reference_quantity);
        else if (From == unit_system::norm_lambda)
            return reference_quantity;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for length from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_length_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(reference_quantity/phys::light_speed);
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0)/reference_quantity;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for area from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_area_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::light_speed *phys::light_speed /
                (reference_quantity*reference_quantity));
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(reference_quantity*reference_quantity);
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for area from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_area_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(reference_quantity*reference_quantity/
                (phys::light_speed*phys::light_speed));
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0/
                (reference_quantity*reference_quantity));
        else
            return static_cast<RealType>(1.0);
    }

        /**
    * This function returns the conversion factor for volume from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_volume_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(
                (phys::light_speed/reference_quantity)*
                (phys::light_speed/reference_quantity)*
                (phys::light_speed/reference_quantity));
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(reference_quantity*reference_quantity*reference_quantity);
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for volume from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_volume_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(
                (reference_quantity/phys::light_speed)*
                (reference_quantity/phys::light_speed)*
                (reference_quantity/phys::light_speed));
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(1.0/
                (reference_quantity*reference_quantity*reference_quantity));
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for time from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_time_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(1.0)/reference_quantity;
        else if (From == unit_system::norm_lambda)
            return reference_quantity*static_cast<RealType>(1.0/phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for time from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_time_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return reference_quantity;
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(phys::light_speed)/reference_quantity;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for rate from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_rate_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return reference_quantity;
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(phys::light_speed)/reference_quantity;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for rate from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_rate_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(1.0)/reference_quantity;
        else if (To == unit_system::norm_lambda)
            return reference_quantity*static_cast<RealType>(1.0/phys::light_speed);
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for electric field from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_E_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::electron_mass*phys::light_speed/
                phys::elementary_charge)*reference_quantity;
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(
                phys::electron_mass*phys::light_speed*phys::light_speed/
                (phys::elementary_charge))/reference_quantity;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for electric field from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_E_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>((phys::elementary_charge)/
                (phys::electron_mass*phys::light_speed))/reference_quantity;
        else if (To == unit_system::norm_lambda)
            return reference_quantity*static_cast<RealType>(
                (phys::elementary_charge)/
                (phys::electron_mass*phys::light_speed*phys::light_speed));
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for magnetic field from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system From, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_B_to_SI_from(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(From == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (From == unit_system::norm_omega)
            return static_cast<RealType>(phys::electron_mass*reference_quantity/
                phys::elementary_charge);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(
                phys::electron_mass*phys::light_speed/
                (reference_quantity*phys::elementary_charge));
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for magnetic field from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity the reference frequency or the ref. wavelength (depending on the unit system) in SI units
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_B_from_SI_to(
        const RealType reference_quantity = static_cast<RealType>(1.0))
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(phys::elementary_charge/
                (phys::electron_mass*reference_quantity));
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(
                (reference_quantity*phys::elementary_charge)/
                (phys::electron_mass*phys::light_speed));
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for energy from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
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

    /**
    * This function returns the conversion factor for energy from SI to
    * any unit system
    *
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_energy_from_SI_to()
    {
        if(To == unit_system::SI)
            return static_cast<RealType>(1.0);
        else if (To == unit_system::norm_omega)
            return static_cast<RealType>(
                1.0/(phys::electron_mass*phys::light_speed*phys::light_speed));
        else if (To == unit_system::norm_lambda)
            return static_cast<RealType>(
                1.0/(phys::electron_mass*phys::light_speed*phys::light_speed));
        else
            return static_cast<RealType>(1.0);
    }

}
}
}
#endif //PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
