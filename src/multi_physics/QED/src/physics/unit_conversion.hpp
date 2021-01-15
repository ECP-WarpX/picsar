#ifndef PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
#define PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION

//Should be included by all the src files of the library
#include "picsar/src/multi_physics/QED/src/qed_commons.h"

//Uses several physical constants
#include "picsar/src/multi_physics/QED/src/physics/phys_constants.h"
//Uses several mathematical constants
#include "picsar/src/multi_physics/QED/src/math/math_constants.h"

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
    * > Natural Heaviside Lorentz units, a unit system frequently used in
    *   particle phyics where c (speed of light) = hbar (reduced Plank constant)
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

    /**
    * All supported unit systems
    */
    enum class unit_system
    {
        SI,
        norm_omega,
        norm_lambda,
        heaviside_lorentz
    };

    /**
    * Reference energy used for the heaviside_lorentz units
    */
    template <typename RealType>
    constexpr RealType heaviside_lorentz_reference_energy = phys::MeV<RealType>;

    /**
    * All possible quantities
    */
    enum class quantity
    {
        mass,
        charge,
        velocity,
        momentum,
        energy,
        length,
        area,
        volume,
        time,
        rate,
        E,
        B
    };

    /**
    * This struct provides a workaround to achieve
    * partial function specialization
    *
    * @tparam Quantity physical quantity to convert
    * @tparam From unit system (enum unit_system) from which to convert
    * @tparam To unit system (enum unit_system) from which to convert
    * @tparam RealType the floating point type of the result
    */
    template<
        quantity Quantity,
        unit_system From,
        unit_system To,
        typename RealType = double
        >
    struct conv
    {
        /**
        * This function returns the conversion factor for any quantity from any unit
        * system to any other unit system
        *
        * @param[in] reference_quantity_from reference quantity needed for From (if needed)
        * @param[in] reference_quantity_to reference quantity needed for To  (if needed)
        * @return the conversion factor
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept;
    };

    //All possible template specializations follow

    // __________________________ mass ____________________________

    template<unit_system UnitSystem, typename RealType>
    struct conv<quantity::mass, UnitSystem, UnitSystem, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/electron_mass<double>);}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/electron_mass<double>);}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            (light_speed<double>*light_speed<double>)/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return electron_mass<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            electron_mass<double>*
            (light_speed<double>*light_speed<double>)/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return electron_mass<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            electron_mass<double>*
            (light_speed<double>*light_speed<double>)/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/
            (light_speed<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/
            (light_speed<double>*light_speed<double>)/
            electron_mass<double>);}
    };

    template<typename RealType>
    struct conv<quantity::mass, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/
            (light_speed<double>*light_speed<double>)/
            electron_mass<double>);}
    };

    // __________________________ charge ____________________________

    template<unit_system UnitSystem, typename RealType>
    struct conv<quantity::charge, UnitSystem, UnitSystem, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/elementary_charge<double>);}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/elementary_charge<double>);}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            sqrt_4_pi_fine_structure<double>/
            elementary_charge<double>);}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return elementary_charge<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return sqrt_4_pi_fine_structure<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return elementary_charge<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return sqrt_4_pi_fine_structure<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                elementary_charge<double>/sqrt_4_pi_fine_structure<double>);}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/sqrt_4_pi_fine_structure<double>);}
    };

    template<typename RealType>
    struct conv<quantity::charge, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/sqrt_4_pi_fine_structure<double>);}
    };

    // __________________________ velocity ____________________________

    template<unit_system UnitSystem, typename RealType>
    struct conv<quantity::velocity, UnitSystem, UnitSystem, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return light_speed<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return light_speed<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return light_speed<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::velocity, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ momentum ____________________________

    template<unit_system UnitSystem, typename RealType>
    struct conv<quantity::momentum, UnitSystem, UnitSystem, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/
            (electron_mass<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/
            (electron_mass<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(light_speed<double>/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                electron_mass<double>*light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(electron_mass<double>*light_speed<double>
            *light_speed<double>/heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                electron_mass<double>*light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(electron_mass<double>*light_speed<double>
            *light_speed<double>/heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            electron_mass<double>/light_speed<double>/light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::momentum, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            electron_mass<double>/light_speed<double>/light_speed<double>);}
    };

    // __________________________ energy ____________________________

    template<unit_system UnitSystem, typename RealType>
    struct conv<quantity::energy, UnitSystem, UnitSystem, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/
            (electron_mass<double>*light_speed<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/
            (electron_mass<double>*light_speed<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(1.0/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                electron_mass<double>*light_speed<double>*light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(electron_mass<double>*light_speed<double>*
            light_speed<double>/heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                electron_mass<double>*light_speed<double>*light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(electron_mass<double>*light_speed<double>
            *light_speed<double>/heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return heaviside_lorentz_reference_energy<double>;}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            electron_mass<double>/light_speed<double>/light_speed<double>);}
    };

    template<typename RealType>
    struct conv<quantity::energy, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            electron_mass<double>/light_speed<double>/light_speed<double>);}
    };

    // __________________________ length ____________________________

    template<typename RealType>
    struct conv<quantity::length, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<RealType>);
            return reference_quantity_to*ff;
        }
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {return math::one<RealType>/reference_quantity_to;}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return light_speed<RealType>/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_to/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return light_speed<RealType>/
            (reference_quantity_from*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/reduced_plank<double>);
            return ff/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<double>);
            return reference_quantity_from*reference_quantity_to*ff;
        }
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_from/reference_quantity_to;}
    };


    template<typename RealType>
    struct conv<quantity::length, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_from*/) noexcept
        {
            constexpr auto ff = static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>));
            return reference_quantity_from*ff;
        }
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>((reduced_plank<double>*light_speed<double>)/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>/
                heaviside_lorentz_reference_energy<double>);
            return ff*reference_quantity_to;
       }
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>*light_speed<double>/
                heaviside_lorentz_reference_energy<double>);
            return ff/reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::length, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ area ____________________________

    template<typename RealType>
    struct conv<quantity::area, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<RealType>);
            return (reference_quantity_to*ff)*(reference_quantity_to*ff);
        }
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {return math::one<RealType>/
            (reference_quantity_to*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>*
            heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>)/
            (reduced_plank<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                light_speed<RealType>*light_speed<RealType>);
            return ff/reference_quantity_from/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return (reference_quantity_to/reference_quantity_from)*
             (reference_quantity_to/reference_quantity_from);}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                light_speed<double>*light_speed<double>);
            return ff/(
                (reference_quantity_from*reference_quantity_to)*
                (reference_quantity_from*reference_quantity_to));}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>/reduced_plank<double>);
            return (ff/reference_quantity_from)*(ff/reference_quantity_from);}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return reference_quantity_from*reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                1.0/light_speed<double>/light_speed<double>);
            return ff*(reference_quantity_from*reference_quantity_to)*
                (reference_quantity_from*reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return (reference_quantity_from/reference_quantity_to)*
                (reference_quantity_from/reference_quantity_to);}
    };


    template<typename RealType>
    struct conv<quantity::area, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>));
            return (reference_quantity_from*ff)*(reference_quantity_from*ff);
        }
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            (reduced_plank<double>*light_speed<double>)*
            (reduced_plank<double>*light_speed<double>)/
            heaviside_lorentz_reference_energy<double>/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>/
                heaviside_lorentz_reference_energy<double>);
            return (ff*reference_quantity_to)*(ff*reference_quantity_to);
       }
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>*light_speed<double>/
                heaviside_lorentz_reference_energy<double>);
            return (ff/reference_quantity_to)*(ff/reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::area, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ volume ____________________________

    template<typename RealType>
    struct conv<quantity::volume, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<double>);
            return (reference_quantity_to*ff)*(reference_quantity_to*ff)*
                (reference_quantity_to*ff);
        }
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {return math::one<RealType>/
            (reference_quantity_to*reference_quantity_to*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>*
            heaviside_lorentz_reference_energy<double>*
            heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>)/
            (reduced_plank<double>*light_speed<double>)/
            (reduced_plank<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                light_speed<RealType>*light_speed<RealType>*light_speed<RealType>);
            return ff/reference_quantity_from/reference_quantity_from/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return (reference_quantity_to/reference_quantity_from)*
             (reference_quantity_to/reference_quantity_from)*
             (reference_quantity_to/reference_quantity_from);}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                light_speed<double>*light_speed<double>*light_speed<double>);
            return ff/(
                (reference_quantity_from*reference_quantity_to)*
                (reference_quantity_from*reference_quantity_to)*
                (reference_quantity_from*reference_quantity_to));}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>/reduced_plank<double>);
            return (ff/reference_quantity_from)*(ff/reference_quantity_from)*
                (ff/reference_quantity_from);}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return reference_quantity_from*reference_quantity_from*reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                1.0/light_speed<double>/light_speed<double>/light_speed<double>);
            return ff*(reference_quantity_from*reference_quantity_to)*
                (reference_quantity_from*reference_quantity_to)*
                (reference_quantity_from*reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return (reference_quantity_from/reference_quantity_to)*
            (reference_quantity_from/reference_quantity_to)*
            (reference_quantity_from/reference_quantity_to);}
    };


    template<typename RealType>
    struct conv<quantity::volume, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>));
            return (reference_quantity_from*ff)*(reference_quantity_from*ff)*
                (reference_quantity_from*ff);
        }
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            (reduced_plank<double>*light_speed<double>)*
            (reduced_plank<double>*light_speed<double>)*
            (reduced_plank<double>*light_speed<double>)/
            heaviside_lorentz_reference_energy<double>/
            heaviside_lorentz_reference_energy<double>/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>/
                heaviside_lorentz_reference_energy<double>);
            return (ff*reference_quantity_to)*(ff*reference_quantity_to)*
                (ff*reference_quantity_to);
       }
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>*light_speed<double>/
                heaviside_lorentz_reference_energy<double>);
            return (ff/reference_quantity_to)*(ff/reference_quantity_to)*
                (ff/reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::volume, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ time ____________________________

    template<typename RealType>
    struct conv<quantity::time, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            return reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {return light_speed<RealType>/reference_quantity_to;}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>/
                reduced_plank<double>);}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_to/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return light_speed<RealType>/
            (reference_quantity_from*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/reduced_plank<double>);
            return ff/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return reference_quantity_from/light_speed<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<double>);
            return reference_quantity_from*reference_quantity_to*ff;
        }
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_from/reference_quantity_to;}
    };


    template<typename RealType>
    struct conv<quantity::time, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(heaviside_lorentz_reference_energy<double>/
            (reduced_plank<double>*light_speed<double>));
            return reference_quantity_from*ff;
        }
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(reduced_plank<double>/
            heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>/
                heaviside_lorentz_reference_energy<double>);
            return ff*reference_quantity_to;
       }
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                reduced_plank<double>*light_speed<double>/
                heaviside_lorentz_reference_energy<double>);
            return ff/reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::time, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ rate ____________________________

    template<typename RealType>
    struct conv<quantity::rate, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {return math::one<RealType>/reference_quantity_to;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_to/light_speed<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(reduced_plank<double>/
                heaviside_lorentz_reference_energy<double>);}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_from/reference_quantity_to;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<double>);
            return ff*(reference_quantity_from*reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>/reduced_plank<double>);
            return reference_quantity_from/ff;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return light_speed<RealType>/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return light_speed<double>/
            (reference_quantity_from*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_to/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                (reduced_plank<double>*light_speed<double>)/
                heaviside_lorentz_reference_energy<double>);
            return ff/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            heaviside_lorentz_reference_energy<double>)/
            reduced_plank<double>;}
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ ,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>/
                reduced_plank<double>);
            return ff/reference_quantity_to;
       }
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>/
                reduced_plank<double>/light_speed<double>);
            return ff*reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::rate, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ E field ____________________________

    template<typename RealType>
    struct conv<quantity::E, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                elementary_charge<double>/electron_mass<double>/
                light_speed<double>);
            return ff/reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                elementary_charge<double>/electron_mass<double>/
                light_speed<double>/light_speed<double>);
            return ff*reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            (elementary_charge<double>*reduced_plank<double>*
            light_speed<double>)/
            heaviside_lorentz_reference_energy<double>/
            heaviside_lorentz_reference_energy<double>/
            sqrt_4_pi_fine_structure<double>);}
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(electron_mass<double>*
                light_speed<double>/elementary_charge<double>);
            return ff*reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_from/reference_quantity_to;}
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<double>);
            return ff*(reference_quantity_from*reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                electron_mass<double>*light_speed<double>*reduced_plank<double>*
                light_speed<double>/heaviside_lorentz_reference_energy<double>/
                heaviside_lorentz_reference_energy<double>/
                sqrt_4_pi_fine_structure<double>);
            return ff*reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                electron_mass<double>*light_speed<double>*light_speed<double>/
                (elementary_charge<double>));
            return ff/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return light_speed<double>/
            (reference_quantity_from*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_to/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(electron_mass<double>*
                light_speed<double>*light_speed<double>*reduced_plank<double>*
                light_speed<double>/
                heaviside_lorentz_reference_energy<double>/
                heaviside_lorentz_reference_energy<double>/
                sqrt_4_pi_fine_structure<double>);
            return ff/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/
                (elementary_charge<double>*reduced_plank<double>*
                light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/reduced_plank<double>/
                light_speed<double>/electron_mass<double>/light_speed<double>);
            return ff /reference_quantity_to;
       }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/reduced_plank<double>/
                light_speed<double>/electron_mass<double>/light_speed<double>/
                light_speed<double>);
            return ff*reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::E, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    // __________________________ B field ____________________________

    template<typename RealType>
    struct conv<quantity::B, unit_system::SI, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::SI, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                elementary_charge<double>/electron_mass<double>);
            return ff/reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::SI, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                elementary_charge<double>/electron_mass<double>/
                light_speed<double>);
            return ff*reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::SI, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
            (elementary_charge<double>*reduced_plank<double>*
            light_speed<double>)*light_speed<double>/
            heaviside_lorentz_reference_energy<double>/
            heaviside_lorentz_reference_energy<double>/
            sqrt_4_pi_fine_structure<double>);}
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_omega, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(electron_mass<double>/
                elementary_charge<double>);
            return ff*reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_omega, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_from/reference_quantity_to;}
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_omega, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(1.0/light_speed<double>);
            return ff*(reference_quantity_from*reference_quantity_to);
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_omega, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                electron_mass<double>*light_speed<double>*reduced_plank<double>*
                light_speed<double>/heaviside_lorentz_reference_energy<double>/
                heaviside_lorentz_reference_energy<double>/
                sqrt_4_pi_fine_structure<double>);
            return ff*reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_lambda, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                electron_mass<double>*light_speed<double>/
                (elementary_charge<double>));
            return ff/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_lambda, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return light_speed<double>/
            (reference_quantity_from*reference_quantity_to);}
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_lambda, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType reference_quantity_to) noexcept
        {return reference_quantity_to/reference_quantity_from;}
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::norm_lambda, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType reference_quantity_from,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {
            constexpr auto ff = static_cast<RealType>(electron_mass<double>*
                light_speed<double>*light_speed<double>*reduced_plank<double>*
                light_speed<double>/
                heaviside_lorentz_reference_energy<double>/
                heaviside_lorentz_reference_energy<double>/
                sqrt_4_pi_fine_structure<double>);
            return ff/reference_quantity_from;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::heaviside_lorentz, unit_system::SI, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/
                (elementary_charge<double>*reduced_plank<double>*
                light_speed<double>*light_speed<double>));}
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::heaviside_lorentz, unit_system::norm_omega, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/reduced_plank<double>/
                light_speed<double>/electron_mass<double>/light_speed<double>);
            return ff /reference_quantity_to;
       }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::heaviside_lorentz, unit_system::norm_lambda, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/,
            const RealType reference_quantity_to) noexcept
        {
            constexpr auto ff = static_cast<RealType>(
                heaviside_lorentz_reference_energy<double>*
                heaviside_lorentz_reference_energy<double>*
                sqrt_4_pi_fine_structure<double>/reduced_plank<double>/
                light_speed<double>/electron_mass<double>/light_speed<double>/
                light_speed<double>);
            return ff*reference_quantity_to;
        }
    };

    template<typename RealType>
    struct conv<quantity::B, unit_system::heaviside_lorentz, unit_system::heaviside_lorentz, RealType>
    {
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        static constexpr RealType fact(
            const RealType /*reference_quantity_from*/ = math::one<RealType>,
            const RealType /*reference_quantity_to*/ = math::one<RealType>) noexcept
        {return math::one<RealType>;}
    };

    //____________________________

    /**
    * Electron rest energy in heaviside_lorentz units
    */
    template <typename RealType>
    constexpr RealType heaviside_lorentz_electron_rest_energy = RealType(
        phys::electron_mass<double>*phys::light_speed<double>*phys::light_speed<double>*
        conv<quantity::energy, unit_system::SI, unit_system::heaviside_lorentz, double>::fact());

    /**
    * Schwinger field in heaviside_lorentz units
    */
    template <typename RealType>
    constexpr RealType heaviside_lorentz_schwinger_field= static_cast<RealType>(
        phys::schwinger_field<double>*conv<quantity::E, unit_system::SI,
            unit_system::heaviside_lorentz, double>::fact());

    /**
    * Elementary charge in heaviside_lorentz units
    */
    template <typename RealType>
    constexpr RealType heaviside_lorentz_elementary_charge = static_cast<RealType>(
        phys::elementary_charge<double>*conv<quantity::charge, unit_system::SI,
            unit_system::heaviside_lorentz, double>::fact());

}
}
}
#endif //PICSAR_MULTIPHYSICS_PHYS_UNIT_CONVERSION
