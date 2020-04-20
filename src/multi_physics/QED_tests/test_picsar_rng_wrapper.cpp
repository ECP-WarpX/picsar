//####### Test module for rng_wrapper ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/rng_wrapper"

//Will automatically define a main for this test
#define BOOST_TEST_DYN_LINK

#include <random>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

#include "rng_wrapper.hpp"

using namespace picsar::multi_physics::utils;

// ------------- Tests --------------

// Warning: random number distributions are not fully checked.
// Only the fact that the extracted number is within the expected
// range is checked.

//***Test constructor & copy

template<typename RealType>
void  rng_stl_wrapper_constructor_copy(uint64_t seed)
{
    std::mt19937_64 rng{seed};

    auto wrp1 = stl_rng_wrapper{seed};
    auto wrp2 = wrp1;

    BOOST_CHECK_EQUAL(
        wrp1.unf<RealType>(0.0,1.0), wrp2.unf<RealType>(0.0,1.0));
    BOOST_CHECK_EQUAL(
        wrp1.exp<RealType>(1.0), wrp2.exp<RealType>(1.0));
    BOOST_CHECK_EQUAL(
        wrp1.poisson<RealType>(1.0), wrp2.poisson<RealType>(1.0));
    BOOST_CHECK_EQUAL(
        wrp1.gaussian<RealType>(0.0,1.0), wrp2.gaussian<RealType>(0.0,1.0));
}

BOOST_AUTO_TEST_CASE( picsar_rng_stl_wrapper_constructor_copy)
{
    rng_stl_wrapper_constructor_copy<double>(2391892344079);
    rng_stl_wrapper_constructor_copy<float>(2391892344079);
}

// *******************************

// ***Test Uniform distribution

template<typename RealType>
void rng_stl_wrapper_unf(uint64_t seed)
{
    auto wrp = stl_rng_wrapper{seed};
    const size_t how_many = 10000;
    const auto a = static_cast<RealType>(-7.0);
    const auto b = static_cast<RealType>(11.1);

    for (size_t d = 0; d <= how_many; ++d){
        const auto qq = wrp.unf(a,b);
        BOOST_TEST( qq >= a);
        BOOST_TEST( qq < b);
    }
}

BOOST_AUTO_TEST_CASE( picsar_rng_stl_wrapper_unf )
{
    rng_stl_wrapper_unf<double>(2391892344079);
    rng_stl_wrapper_unf<float>(2391892344079);
}

// *******************************

// ***Test exponential distribution

template<typename RealType>
void rng_stl_wrapper_exp(uint64_t seed)
{
    auto wrp = stl_rng_wrapper{seed};
    const size_t how_many = 10000;
    const auto l = static_cast<RealType>(1.0);

    for (size_t d = 0; d <= how_many; ++d){
        const auto qq = wrp.exp(l);
        BOOST_TEST( qq >= 0.0);
    }
}

BOOST_AUTO_TEST_CASE( picsar_rng_stl_wrapper_exp )
{
    rng_stl_wrapper_exp<float>(2391892344079);
    rng_stl_wrapper_exp<double>(2391892344079);
}

// *******************************

// ***Test Poisson distribution

template<typename RealType>
void rng_stl_wrapper_poisson(uint64_t seed)
{
    auto wrp = stl_rng_wrapper{seed};
    const size_t how_many = 10000;
    const auto l = static_cast<RealType>(1.0);

    for (size_t d = 0; d <= how_many; ++d){
        const auto qq = wrp.poisson(l);
        BOOST_TEST( qq >= 0);
    }
}

BOOST_AUTO_TEST_CASE( picsar_rng_stl_wrapper_poisson )
{
    rng_stl_wrapper_poisson<float>(2391892344079);
    rng_stl_wrapper_poisson<double>(2391892344079);
}

// *******************************

// ***Test gaussian distribution

template<typename RealType>
void rng_stl_wrapper_gaussian(uint64_t seed)
{
    auto wrp = stl_rng_wrapper{seed};
    const size_t how_many = 10000;
    const auto l = static_cast<RealType>(1000.0);
    const auto sigma = static_cast<RealType>(0.1);

    for (size_t d = 0; d <= how_many; ++d){
        const auto qq = wrp.gaussian(l, sigma);
        BOOST_TEST( qq > 0.0); //qq < 0 is VERY unlikely
    }
}

BOOST_AUTO_TEST_CASE( picsar_rng_stl_wrapper_gaussian )
{
    rng_stl_wrapper_gaussian<float>(2391892344079);
    rng_stl_wrapper_gaussian<double>(2391892344079);
}

// *******************************
