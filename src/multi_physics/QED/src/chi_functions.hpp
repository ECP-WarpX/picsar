#ifndef __PICSAR_MULTIPHYSICS_CHI_FUNCTIONS__
#define __PICSAR_MULTIPHYSICS_CHI_FUNCTIONS__

//This .hpp file contais the functions to calcuate the chi parameter for
//photons and leptons

#include <cmath>

//Should be included by all the src files of the library
#include "qed_commons.h"

//Uses vectors
#include "vec_functions.hpp"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //chi for photons
        //lambda will be ignored if compiled with SI units
        template<typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL calc_chi_photon(
            _REAL px, _REAL py, _REAL pz,
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL lambda = static_cast<_REAL>(1.0)
        );

        //chi for leptons
        //lambda will be ignored if compiled with SI units
        template<typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL calc_chi_leptons(
            _REAL px, _REAL py, _REAL pz,
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL lambda = static_cast<_REAL>(1.0)
        );
    }
}

//############################################### Implementation

//chi for photons
//NB: constants used in this function (e.g. __c instead of light_speed)
//take into account the unit convention
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::calc_chi_photon(
    _REAL px, _REAL py, _REAL pz,
    _REAL ex, _REAL ey, _REAL ez,
    _REAL bx, _REAL by, _REAL bz,
    _REAL lambda
)
{
#ifdef PXRMP_WITH_SI_UNITS //Enforces lambda=1.0 if SI units are used
    lambda = static_cast<_REAL>(1.0);
#endif

    vec3<_REAL> p{px, py, pz};
    vec3<_REAL> em_e{ex, ey, ez};
    vec3<_REAL> em_b{bx, by, bz};

    vec3<_REAL> p_unit = p / norm(p);
    vec3<_REAL> em_eperp = em_e - dot(p_unit,em_e)*p_unit;
    _REAL mod = norm(em_eperp + cross(p_unit*static_cast<_REAL>(__c), em_b));

    return mod*p/(static_cast<_REAL>(__emass*__c*__schwinger)*lambda);
}

//chi for leptons
//NB: constants used in this function (e.g. __c instead of light_speed)
//take into account the unit convention
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::calc_chi_leptons(
    _REAL px, _REAL py, _REAL pz,
    _REAL ex, _REAL ey, _REAL ez,
    _REAL bx, _REAL by, _REAL bz,
    _REAL lambda
)
{
#ifdef PXRMP_WITH_SI_UNITS //Enforces lambda=1.0 if SI units are used
    lambda = static_cast<_REAL>(1.0);
#endif

    return 0.0;
}



#endif //__PICSAR_MULTIPHYSICS_CHI_FUNCTIONS__
