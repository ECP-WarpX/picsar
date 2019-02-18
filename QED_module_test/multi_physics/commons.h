#ifndef __PMP_COMMONS__
#define __PMP_COMMONS__

#include<cstdint>

//Common definitions for the multi_physics library (header only)

namespace picsar{
    namespace multi_physics{
        enum particle_type(electrons, positrons, photons, ions);

        const double _km = 1.0e3;
        const double _m = 1.0;
        const double _cm = 1.0e-2;
        const double _mm = 1.0e-3;
        const double _um = 1.0e-6;
        const double _nm = 1.0e-9;

        const double classical_electron_radius = 2.8179403267e-13 * _cm;

        //some checks on the size of the fundamental types
        static_assert(sizeof(int) == 4, "int type is not 32 bit!");
        static_assert(sizeof(double) == 8, "double type is not 64 bit!");
    }
}

#endif
