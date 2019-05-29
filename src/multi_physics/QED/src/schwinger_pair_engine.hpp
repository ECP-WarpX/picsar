#ifndef __PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE__
#define __PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE__

//This .hpp file contais the implementation of the
//Schwinger pair engine (as described in Gonoskov et al. PRE 92, 023305 2015
// and Banerjee et al. PRA 98, 032121 2018)

//Uses pairs
#include <tuple>

//Should be included by all the src files of the library
#include "qed_commons.h"

//Uses picsar arrays
#include "picsar_array.hpp"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{
        template<typename _REAL, class _RNDWRAP>
        class schwinger_pair_engine
        {
        public:
            //A random number generatator has to be passed by move.
            //The RNG can be ANY object implementing the functions
            //_REAL unf (_REAL a, _REAL b)
            //and
            //_REAL exp (_REAL l)
            //The constructor can accept a lambda parameter.
            //It is ignored if the SI units option is selected
            schwinger_pair_engine(_RNDWRAP&& rng,
            _REAL lambda = static_cast<_REAL>(1.0));

            //Copy constructor
            schwinger_pair_engine(const schwinger_pair_engine& other);

            //Move constructor
            schwinger_pair_engine(schwinger_pair_engine&& other);

            //Getter & setter for lambda
            _REAL get_lambda() const;
            void set_lambda(_REAL lambda);

            //This function determines if a pair has been generated in a given
            //cell and its statistical weight.
            //It returns a bool (true if a pair has been generated)
            //and a _REAL (the statistical weight of the new pair)
            //It requires to provide the fields, the cell size and dt
            //Use this function if the probability to generate a pair is << 1
            PXRMP_FORCE_INLINE
            std::pair<bool, _REAL> generate_pairs_single(
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL dx, _REAL dy, _REAL dz,
            _REAL dt
            );

            //______________________GPU
            //Same as above, but with GPU use directly this one!
            PXRMP_GPU
            PXRMP_FORCE_INLINE
            static
            void
            internal_generate_pairs_single(
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL dx, _REAL dy, _REAL dz,
            _REAL dt,
            bool* has_event_happend, _REAL* weight,
            _REAL lambda, _REAL unf_zero_one_minus_epsi
            );


            //This function determines how many pairs have been generated in a given
            //cell and their statistical weight.
            //It returns a size_t (the number of the generated particles)
            //and a _REAL (the statistical weight of the new pair)
            //It requires to provide the fields, the cell size and dt
            //Use this function if the probability to generate a pair is large
            PXRMP_FORCE_INLINE
            std::pair<size_t, _REAL> generate_pairs_multiple(
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL dx, _REAL dy, _REAL dz,
            _REAL dt
            );

            //______________________GPU
            //Same as above, but with GPU use directly this one!
            PXRMP_GPU
            PXRMP_FORCE_INLINE
            static
            void
            internal_generate_pairs_multiple(
            _REAL ex, _REAL ey, _REAL ez,
            _REAL bx, _REAL by, _REAL bz,
            _REAL dx, _REAL dy, _REAL dz,
            _REAL dt,
            size_t* how_many, _REAL* weight,
            _REAL lambda, _REAL unf_zero_one_minus_epsi
            );

            //This function computes 3 random numbers between 0 and 1.
            //They should be used to initialize the position of the pair.
            PXRMP_FORCE_INLINE
            picsar_array<_REAL, 3>
            get_new_pair_position();

        private:

            _REAL lambda;
            //The only requrement for the RNG is to be able to provide unf(a,b) and
            //exp(l)
            _RNDWRAP rng;

            //Some handy constants
            static constexpr _REAL zero = static_cast<_REAL>(0.0);
            static constexpr _REAL one = static_cast<_REAL>(1.0);


        };
    }
}

//###############################################

template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
schwinger_pair_engine
(_RNDWRAP&& rng, _REAL lambda):
    lambda{lambda}, rng{std::move(rng)}
{
    //This enforces lambda=1 if SI units are used.
#ifdef PXRMP_WITH_SI_UNITS
    lambda = static_cast<_REAL>(1.0);
#endif
}

//Copy constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
schwinger_pair_engine(const schwinger_pair_engine& other):
    lambda(other.lambda), rng(other.rng)
    {}

//Move constructor
template<typename _REAL, class _RNDWRAP>
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
schwinger_pair_engine(schwinger_pair_engine&& other):
    lambda(std::move(other.lambda)), rng(std::move(other.rng))
    {}


//Getter for lambda
template<typename _REAL, class _RNDWRAP>
_REAL picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
get_lambda() const
{
    return lambda;
}

//Setter for lambda
template<typename _REAL, class _RNDWRAP>
void picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
set_lambda
#ifdef PXRMP_WITH_NORMALIZED_UNITS
(_REAL lambda)
{
    this->lambda = lambda;
}
#else
(_REAL){} //Do nothing
#endif


//This function determines if a pair has been generated in a given
//cell and its statistical weight.
//It returns a bool (true if a pair has been generated)
//and a _REAL (the statistical weight of the new pair)
//It requires to provide the fields, the cell size and dt
//Use this function if the probability to generate a pair is << 1
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
std::pair<bool, _REAL>
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
generate_pairs_single(
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dx, _REAL dy, _REAL dz,
_REAL dt
)
{
    std::pair<bool, _REAL> res{false, zero};

    internal_generate_pairs_single(ex, ey, ez, bx, by, bz, dt,
        &res.first, &res.second, lambda, rng.unf(zero, one));

    return res;
}


//______________________GPU
//Same as above, but with GPU use directly this one!
//Returns false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
void
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
internal_generate_pairs_single(
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dx, _REAL dy, _REAL dz,
_REAL dt,
bool* has_event_happend, _REAL* weight,
_REAL lambda, _REAL unf_zero_one_minus_epsi
)
{
    //All the magic happens here!
}

//This function determines how many pairs have been generated in a given
//cell and their statistical weight.
//It returns a size_t (the number of the generated particles)
//and a _REAL (the statistical weight of the new pair)
//It requires to provide the fields, the cell size and dt
//Use this function if the probability to generate a pair is large
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
std::pair<size_t, _REAL>
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
generate_pairs_multiple(
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dx, _REAL dy, _REAL dz,
_REAL dt
)
{
    std::pair<size_t, _REAL> res{0, zero};

    internal_generate_pairs_multiple(ex, ey, ez, bx, by, bz, dt,
        &res.first, &res.second, lambda, rng.unf(zero, one));

    return res;
}

//______________________GPU
//Same as above, but with GPU use directly this one!
//Returns false if errors occur
template<typename _REAL, class _RNDWRAP>
PXRMP_GPU
PXRMP_FORCE_INLINE
void
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
internal_generate_pairs_multiple(
_REAL ex, _REAL ey, _REAL ez,
_REAL bx, _REAL by, _REAL bz,
_REAL dx, _REAL dy, _REAL dz,
_REAL dt,
size_t* how_many, _REAL* weight,
_REAL lambda, _REAL unf_zero_one_minus_epsi
)
{
    //All the magic happens here!
}


//This function computes 3 random numbers between 0 and 1.
//They should be used to initialize the position of the pair.
template<typename _REAL, class _RNDWRAP>
PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_array<_REAL, 3>
picsar::multi_physics::schwinger_pair_engine<_REAL, _RNDWRAP>::
get_new_pair_position()
{
    picsar_array<_REAL, 3> res;
    for(auto& val: res)
        val = rng.unf(zero, one);
    return res;
}


#endif //__PICSAR_MULTIPHYSICS_SCHWINGER_PAIR_ENGINE__
