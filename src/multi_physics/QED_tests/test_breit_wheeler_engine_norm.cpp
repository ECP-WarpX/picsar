//####### Test module for Breit Wheeler engine (SI units) ######################

//Define Module name
 #define BOOST_TEST_MODULE "Nonliner Breit Wheeler engine (SI)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//SI units are used for this test
#define PXRMP_USE_NORMALIZED_UNITS
#include "breit_wheeler_engine.hpp"

using namespace picsar::multi_physics;

//Helper function
template<class RNG, typename REAL, typename WHATEVER>
breit_wheeler_engine<REAL, RNG> get_bw_set_lambda(int64_t seed, WHATEVER lambda)
{
    stl_rng_wrapper rng{seed};
    auto bw_engine =  breit_wheeler_engine<REAL, stl_rng_wrapper>{std::move(rng)};
    bw_engine.set_lambda(static_cast<REAL>(lambda));
    return bw_engine;
}

// ------------- Tests --------------

//Test get/set lambda for breit_wheeler_engine (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_double_1 )
{
    auto bw_engine = get_bw_set_lambda<stl_rng_wrapper, double>
        (390109317, 800.0*si_nanometer);
    BOOST_CHECK_EQUAL( 800.0*si_nanometer, bw_engine.get_lambda());
}

//Test get/set lambda for breit_wheeler_engine (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_single_1 )
{
    auto bw_engine = get_bw_set_lambda<stl_rng_wrapper, float>
        (390109317, 800.0*si_nanometer);
    BOOST_CHECK_EQUAL( static_cast<float>(800.0*si_nanometer),
        bw_engine.get_lambda());
}

// ------------- optical depth --------------

//Test get new optical depth (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_opt_double_1 )
{
    auto bw_engine = get_bw_set_lambda<stl_rng_wrapper, double>
        (390109317, 800.0*si_nanometer);
    BOOST_TEST ( bw_engine.get_optical_depth() >= 0.0 );
}

//Test get new optical depth (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_opt_single_1 )
{
    auto bw_engine = get_bw_set_lambda<stl_rng_wrapper, float>
        (390109317, 800.0*si_nanometer);
    BOOST_TEST( bw_engine.get_optical_depth() >= 0.0 );
}
