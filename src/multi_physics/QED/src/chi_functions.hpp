#ifndef __PICSAR_MULTIPHYSICS_CHI_FUNCTIONS__
#define __PICSAR_MULTIPHYSICS_CHI_FUNCTIONS__

//This .hpp file contains the functions to calcuate the chi parameter for
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
        _REAL chi_photon(
            _REAL px, _REAL py, _REAL pz,
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL lambda = static_cast<_REAL>(1.0)
        );

        //chi for leptons
        //lambda will be ignored if compiled with SI units
        template<typename _REAL>
        PXRMP_FORCE_INLINE
        _REAL chi_lepton(
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
_REAL picsar::multi_physics::chi_photon(
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

    _REAL norm_p = norm(p);
    if(norm_p == static_cast<_REAL>(0.0))
        return static_cast<_REAL>(0.0);

    vec3<_REAL> p_unit = p / norm_p;
    vec3<_REAL> em_eperp = em_e - dot(p_unit,em_e)*p_unit;
    _REAL mod = norm(em_eperp + cross(p_unit*static_cast<_REAL>(__c), em_b));

    _REAL coeff = static_cast<_REAL>(1.0)/
        (static_cast<_REAL>(__emass*__c*__schwinger)*lambda);

    return mod*norm_p*coeff;
}

//chi for leptons
//NB: constants used in this function (e.g. __c instead of light_speed)
//take into account the unit convention
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL picsar::multi_physics::chi_lepton(
    _REAL px, _REAL py, _REAL pz,
    _REAL ex, _REAL ey, _REAL ez,
    _REAL bx, _REAL by, _REAL bz,
    _REAL lambda
)
{
#ifdef PXRMP_WITH_SI_UNITS //Enforces lambda=1.0 if SI units are used
    lambda = static_cast<_REAL>(1.0);
#endif

    _REAL one = static_cast<_REAL>(1.0);

    vec3<_REAL> p{px, py, pz};
    vec3<_REAL> em_e{ex, ey, ez};
    vec3<_REAL> em_b{bx, by, bz};

    _REAL norm_p2 = norm2(p);
    if(norm_p2 == static_cast<_REAL>(0.0))
        return static_cast<_REAL>(0.0);

    _REAL norm_p = sqrt(norm_p2);
    vec3<_REAL> p_unit = p / norm_p;

    //For gamma_2, writing the operations like this is better for single
    //precision (found with some tests).
    _REAL gamma_2 = one +
        (norm_p/static_cast<_REAL>(__emass*__c))*
        (norm_p/static_cast<_REAL>(__emass*__c));
    _REAL gamma = sqrt(gamma_2);

    _REAL beta = sqrt(one-one/gamma_2);
    vec3<_REAL> beta_vec = beta * p_unit;

    _REAL beta_dot_e = dot(beta_vec, em_e);
    _REAL beta_dot_e_2 = beta_dot_e*beta_dot_e;
    _REAL e_plus_v_cross_b_2 = norm2(em_e +
        cross(beta_vec * static_cast<_REAL>(__c), em_b));

    _REAL coeff = static_cast<_REAL>(1.0)/
        (static_cast<_REAL>(__schwinger)*lambda);

    return gamma*sqrt(fabs(beta_dot_e_2-e_plus_v_cross_b_2))*coeff;

}



#endif //__PICSAR_MULTIPHYSICS_CHI_FUNCTIONS__
