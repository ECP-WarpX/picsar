#ifndef PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
#define PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION

#include "../qed_commons.h"
#include "phys_constants.h"
#include "../math/math_constants.h"

#include <cmath>

namespace picsar{
namespace multi_physics{
namespace phys{

    /**
    * The interface of the PICSAR QED library supports 4 unit systems
    * (Internally all the calculations are performed in Natural
    * Heaviside-Lorentz units):
    *
    * > SI units (International System of Units)
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
    * > Natural Heaviside Lorentz units, a unit system frequently used in particle phyics
    *   where c (speed of light) = hbar (reduced Plank constant) =
    *   = epsilon_0 (vacuum permittivity) = mu_0 (vacuum_permeability) = 1.
    *   In addition, we choose to measure the energy in MeV units
    *   (1 MeV = 1.0e6 * e0 * 1 V (Volt) = 1.60218e-13 J ).
    *   As a consequence:
    *   - mass is in MeV units (m --> mc^2 [MeV])
    *   - length is in 1/MeV units (l --> l / (hbar*c) [MeV^-1])
    *   - time is in 1/MeV units (t --> t / (hbar) [MeV^-1])
    *   - momentum is in MeV units (p --> p*c [MeV])
    *   - energy is simply expressed in MeV
    *   - velocity is normalized with respect to the speed of light
    *   - electric charge is dimensionless (q --> q*sqrt(4.0*pi*fine_structure)/e0)
    *   - electric field is in MeV^2 units (E --> e0*hbar*c*E/sqrt(4*pi*fine_structure) [MeV^2])
    *   - magnetic field is in MeV^2 (B --> e0*hbar*c^2*E/sqrt(4*pi*fine_structure) [MeV^2])
    *
    * The helper functions provided by this module are used internally by the
    * library to convert the input values of a function to Natural units and to
    * convert the output from Natural units to the desired unit system.
    *
    * Inside these helper functions, calculations of constexpr constants are
    * carried out in double precision.
    */
    enum unit_system
    {
        SI,
        norm_omega,
        norm_lambda,
        heaviside_lorentz
    };

    template <typename RealType>
    constexpr RealType heaviside_lorentz_reference_energy = phys::MeV<RealType>;

    /**
    * This function returns the conversion factor for mass from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system From, typename RealType = double>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_mass_to_SI_from() noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return electron_mass<RealType>;
        else if (From == unit_system::norm_lambda)
            return electron_mass<RealType>;
        else if (From == unit_system::heaviside_lorentz)
            return static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>/
                (light_speed<double>*light_speed<double>));
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for mass from any unit
    * system to any other unit system.
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType = double>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_mass_from_to()noexcept
    {
            return fact_mass_to_SI_from<From,RealType>()/
            fact_mass_to_SI_from<To,RealType>();
    }

    /**
    * This function returns the conversion factor for charge from any unit
    * system to SI
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system From, typename RealType = double>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_charge_to_SI_from() noexcept
    {
        using namespace phys;
        using namespace math;

        if (From == unit_system::norm_omega)
            return elementary_charge<RealType>;
        else if (From == unit_system::norm_lambda)
            return elementary_charge<RealType>;
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res = static_cast<RealType>(
                elementary_charge<double>/sqrt_4_pi_fine_structure<double>);
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for charge from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_charge_from_to() noexcept
    {
        return fact_charge_to_SI_from<From,RealType>()/
            fact_charge_to_SI_from<To,RealType>();
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
    constexpr RealType fact_velocity_to_SI_from() noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return light_speed<RealType>;
        else if (From == unit_system::norm_lambda)
            return light_speed<RealType>;
        else if (From == unit_system::heaviside_lorentz)
            return light_speed<RealType>;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for velocity from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_velocity_from_to() noexcept
    {
        return fact_velocity_to_SI_from<From,RealType>()/
            fact_velocity_to_SI_from<To,RealType>();
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
    constexpr RealType fact_momentum_to_SI_from() noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return static_cast<RealType>(
                electron_mass<double>*light_speed<double>);
        else if (From == unit_system::norm_lambda)
            return static_cast<RealType>(
                electron_mass<double>*light_speed<double>);
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res =
                static_cast<RealType>(
                    heaviside_lorentz_reference_energy<double>/light_speed<double>);
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for momentum from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_momentum_from_to() noexcept
    {
        return fact_momentum_to_SI_from<From,RealType>()/
            fact_momentum_to_SI_from<To,RealType>();
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
        using namespace phys;
        if (From == unit_system::norm_omega){
            constexpr auto res = static_cast<RealType>(
                electron_mass<double>*light_speed<double>*light_speed<double>);
            return res;
        }
        else if (From == unit_system::norm_lambda)
        {
            constexpr auto res = static_cast<RealType>(
                electron_mass<double>*light_speed<double>*light_speed<double>);
            return res;
        }
        else if (From == unit_system::heaviside_lorentz)
            return heaviside_lorentz_reference_energy<RealType>;
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for energy from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_energy_from_to() noexcept
    {
        return fact_energy_to_SI_from<From,RealType>()/
            fact_energy_to_SI_from<To,RealType>();
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return light_speed<RealType>/reference_quantity;
        else if (From == unit_system::norm_lambda)
            return reference_quantity;
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res = static_cast<RealType>(
                reduced_plank<double>*light_speed<double>/
                heaviside_lorentz_reference_energy<double>);
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for length from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_length_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_length_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_length_to_SI_from<To,RealType>(reference_quantity_to);
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return (light_speed<RealType>/reference_quantity)*
                (light_speed<RealType>/reference_quantity);
        else if (From == unit_system::norm_lambda)
            return reference_quantity*reference_quantity;
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res = static_cast<RealType>(
                (reduced_plank<double>*light_speed<double>/
                    heaviside_lorentz_reference_energy<double>)*
                (reduced_plank<double>*light_speed<double>/
                    heaviside_lorentz_reference_energy<double>));
                return res;
            }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for area from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_area_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_area_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_area_to_SI_from<To,RealType>(reference_quantity_to);
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return (light_speed<RealType>/reference_quantity)*
                (light_speed<RealType>/reference_quantity)*
                (light_speed<RealType>/reference_quantity);
        else if (From == unit_system::norm_lambda)
            return reference_quantity*
                reference_quantity*
                reference_quantity;
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res = static_cast<RealType>(
                (reduced_plank<double>*light_speed<double>/
                    heaviside_lorentz_reference_energy<double>)*
                (reduced_plank<double>*light_speed<double>/
                    heaviside_lorentz_reference_energy<double>)*
                (reduced_plank<double>*light_speed<double>/
                    heaviside_lorentz_reference_energy<double>));
                return res;
            }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for volume from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_volume_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_volume_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_volume_to_SI_from<To,RealType>(reference_quantity_to);
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return static_cast<RealType>(1.0)/reference_quantity;
        else if (From == unit_system::norm_lambda)
            return reference_quantity/light_speed<RealType>;
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res =
                static_cast<RealType>(reduced_plank<double>/
                    heaviside_lorentz_reference_energy<double>);
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for time from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_time_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_time_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_time_to_SI_from<To,RealType>(reference_quantity_to);
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;

        if (From == unit_system::norm_omega)
            return reference_quantity;
        else if (From == unit_system::norm_lambda)
            return light_speed<RealType>/reference_quantity;
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res =
                static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
                    reduced_plank<double>);
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for rate from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_rate_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_rate_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_rate_to_SI_from<To,RealType>(reference_quantity_to);
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;
        using namespace math;

        if (From == unit_system::norm_omega){
            constexpr auto fact = static_cast<RealType>(
                electron_mass<double>*light_speed<double>/elementary_charge<double>);
            return fact*reference_quantity;
        }
        else if (From == unit_system::norm_lambda){
            constexpr auto fact = static_cast<RealType>(
                electron_mass<double>*light_speed<double>*light_speed<double>/
                (elementary_charge<double>));
            return fact/reference_quantity;
        }
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/
                (elementary_charge<double>*reduced_plank<double>*
                light_speed<double>));
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for electric field from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_E_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_E_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_E_to_SI_from<To,RealType>(reference_quantity_to);
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
        const RealType reference_quantity = static_cast<RealType>(1.0)) noexcept
    {
        using namespace phys;
        using namespace math;

        if (From == unit_system::norm_omega){
            constexpr auto fact = static_cast<RealType>(
                    electron_mass<double>/elementary_charge<double>);
            return fact * reference_quantity;
        }
        else if (From == unit_system::norm_lambda){
            constexpr auto fact =
                static_cast<RealType>(electron_mass<double>*light_speed<double>/
                    elementary_charge<double>);
            return fact / reference_quantity;
        }
        else if (From == unit_system::heaviside_lorentz){
            constexpr auto res = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/
                (elementary_charge<double>*reduced_plank<double>*
                light_speed<double>*light_speed<double>));
            return res;
        }
        else
            return static_cast<RealType>(1.0);
    }

    /**
    * This function returns the conversion factor for magnetic field from any unit
    * system to any other unit system
    *
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) to which to convert
    * @tparam RealType the floating point type of the result
    * @param reference_quantity_from the reference frequency or the ref. wavelength for From in SI units
    * @param reference_quantity_to the reference frequency or the ref. wavelength for To in SI units
    * @return the conversion factor
    */
    template<unit_system From, unit_system To, typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    constexpr RealType fact_B_from_to(
        const RealType reference_quantity_from = static_cast<RealType>(1.0),
        const RealType reference_quantity_to = static_cast<RealType>(1.0)) noexcept
    {
        return fact_B_to_SI_from<From,RealType>(reference_quantity_from)/
            fact_B_to_SI_from<To,RealType>(reference_quantity_to);
    }
}
}
}
#endif //PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
