#ifndef __PICSAR_MULTIPHYSICS_QED_COMMONS__
#define __PICSAR_MULTIPHYSICS_QED_COMMONS__

//This header file contains usefuls definitions (macros, physical constants...)
//Should be included by all the source files of the QED library

//############################################## Compiler specific macros ######

        //Restrict qualifier (compiler specific!)
        #ifdef _WIN32
            #define PXRMP_RESTRICT __restrict
        #else
            #define PXRMP_RESTRICT __restrict__
        #endif

        //Force inline pragmas (compiler specific!)
        #if defined(__CUDA_ARCH__)
            #define PXRMP_FORCE_INLINE __forceinline__
        #elif defined(__INTEL_COMPILER)
            #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__clang__)
            #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__GNUC__)
            #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__ibmxl__)
            #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #else
            #define PXRMP_FORCE_INLINE inline
        #endif

//##############################################################################

namespace picsar{
    namespace multi_physics{

//######################## Physical & mathematical constants (S.I.) ############

        //Double precision
        const double pi = 3.141592653589793;

        const double electron_mass = 9.10938356e-31;
        const double elementary_charge = 1.6021766208e-19;
        const double light_speed = 299792458;
        const double reduced_plank = 1.054571800e-34;
        const double vacuum_permittivity =  8.854187817e-12;
        const double fine_structure =  0.0072973525664;

        const double classical_electron_radius =
        elementary_charge*elementary_charge /
        (4.0*pi*vacuum_permittivity*electron_mass*light_speed*light_speed);

        //Single precision
        const float flt_pi = static_cast<float>(pi);

        const float flt_electron_mass =
            static_cast<float>(electron_mass);
        const float flt_elementary_charge =
            static_cast<float>(elementary_charge);
        const float flt_light_speed =
            static_cast<float>(light_speed);
        const float flt_reduced_plank =
            static_cast<float>(flt_reduced_plank);
        const float flt_vacuum_permittivity =
            static_cast<float>(flt_vacuum_permittivity);
        const float flt_fine_structure =
            static_cast<float>(flt_fine_structure);

        const float flt_classical_electron_radius =
        static_cast<float>(flt_classical_electron_radius);


//##############################################################################

    }
}

#endif// __PICSAR_MULTIPHYSICS_QED_COMMONS__
