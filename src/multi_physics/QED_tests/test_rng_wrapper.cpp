//####### Test module for rng_wrapper ####################################

//Define Module name
 #define BOOST_TEST_MODULE "RNG wrapper"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

#include <cstdint>
#include <random>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS
#include "rng_wrapper.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//***STL***

//***Constructors

//Test STL rng_wrapper constructors generic
template<typename T>
void  rng_stl_wrapper_constructors(uint64_t seed)
{
    std::mt19937_64 rng{seed};

    stl_rng_wrapper<T> wrp1{seed};
    stl_rng_wrapper<T> wrp2(move(rng));

    BOOST_CHECK_EQUAL( wrp1.unf(0.0,1.0), wrp2.unf(0.0,1.0));
}

//Test STL rng_wrapper constructors(double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_constructors_double )
{
    rng_stl_wrapper_constructors<double>(2391892344079);
}

//Test STL rng_wrapper constructors(single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_constructors_single )
{
    rng_stl_wrapper_constructors<float>(2391892344079);
}

//***Uniform distribution

//Test STL rng_wrapper unf generic
template<typename T>
void  rng_stl_wrapper_unf(uint64_t seed)
{
    stl_rng_wrapper<T> wrp{seed};
    size_t how_many = 10000;
    T a = static_cast<T>(-7.0);
    T b = static_cast<T>(11.1);

    for (size_t d = 0; d <= how_many; ++d){
        T qq = wrp.unf(a,b);
        BOOST_TEST( qq >= a);
        BOOST_TEST( qq < b);
    }
}

//Test STL rng_wrapper unf (double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_unf_double )
{
    rng_stl_wrapper_unf<double>(2391892344079);
}

//Test STL rng_wrapper unf (single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_unf_single )
{
    rng_stl_wrapper_unf<float>(2391892344079);
}


//Test STL rng_wrapper exp generic
template<typename T>
void  rng_stl_wrapper_exp(uint64_t seed)
{
    stl_rng_wrapper<T> wrp{seed};
    size_t how_many = 10000;
    T l = static_cast<T>(1.0);

    for (size_t d = 0; d <= how_many; ++d){
        T qq = wrp.exp(l);
        BOOST_TEST( qq >= 0.0);
    }
}


//Test STL rng_wrapper exp (double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_exp_double )
{
    rng_stl_wrapper_exp<double>(2391892344079);
}

//Test STL rng_wrapper exp (single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_exp_single )
{
    rng_stl_wrapper_exp<float>(2391892344079);
}

//***Exponential distribution
