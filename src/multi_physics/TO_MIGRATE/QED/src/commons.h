#ifndef __PMP_COMMONS__
#define __PMP_COMMONS__

#include<cstdint>
#include<cmath>
#include<string>
#include<fstream>
#include<iostream>

//Common definitions for the multi_physics library (header only)

namespace picsar{
    namespace multi_physics{

        //Restrict qualifier (compiler specific!)
        #ifdef _WIN32
            #define PXR_RESTRICT __restrict
        #else
            #define PXR_RESTRICT __restrict__
        #endif

        //Force inline pragmas (compiler specific!)
        #if defined(__CUDA_ARCH__)
            #define PXR_FORCE_INLINE __forceinline__
        #elif defined(__INTEL_COMPILER)
            #define PXR_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__clang__)
            #define PXR_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__GNUC__)
            #define PXR_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__ibmxl__)
            #define PXR_FORCE_INLINE inline __attribute__((always_inline))
        #else
            #define PXR_FORCE_INLINE inline
        #endif


        //Parameters in S.I. units
        const double electron_mass = 9.10938356e-31;
        const double elementary_charge = 1.6021766208e-19;
        const double speed_of_light = 299792458;
        const double reduced_plank = 1.054571800e-34;
        const double vacuum_permittivity =  8.854187817e-12;
        const double fine_structure_constant =  0.0072973525664;

        const double BIG_POSITIVE_DOUBLE = 1.0e30;

        const double classical_electron_radius = elementary_charge*elementary_charge /
                    (4.0*M_PI*vacuum_permittivity*electron_mass*speed_of_light*speed_of_light);

        //Calculate Schwinger field in CODE units given lambda in SI units
        inline double calc_schwinger_given_lambda(double lambda){
             return electron_mass*speed_of_light*lambda / (reduced_plank * 2.0 * M_PI);
        }

        const double _Gm = 1.0e9;
        const double _Mm = 1.0e6;
        const double _km = 1.0e3;
        const double _m = 1.0;
        const double _cm = 1.0e-2;
        const double _mm = 1.0e-3;
        const double _um = 1.0e-6;
        const double _nm = 1.0e-9;
        const double _pm = 1.0e-12;
        const double _fm = 1.0e-15;

        enum particle_type {electrons, positrons, photons, ions};

        //some checks on the size of the fundamental types
        static_assert(sizeof(int) == 4, "int type is not 32 bit!");
        static_assert(sizeof(double) == 8, "double type is not 64 bit!");

        //Message function
        inline void message(const std::string& msg, std::ostream* diag){
            if(diag != nullptr)
                *diag << msg << "\n";
        }

        //Error Message function
        inline void err(const std::string& msg){
            std::cerr << msg << "\n";
        }


    }
}

#endif
