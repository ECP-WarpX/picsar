//####### Test module for Breit Wheeler engine (SI units) ######################

//Define Module name
 #define BOOST_TEST_MODULE "Nonliner Breit Wheeler engine (SI)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include <memory>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//SI units are used for this test
#define PXRMP_USE_SI_UNITS
#include "breit_wheeler_engine.hpp"

#include "rng_wrapper.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Test get/set lambda for breit_wheeler_engine (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_double_1 )
{
    stl_rng_wrapper rng{390109317};

    double lambda = 800.0*si_nanometer;
    breit_wheeler_engine<double, stl_rng_wrapper> bw_engine(std::move(rng));

    bw_engine.set_lambda(lambda);
    BOOST_CHECK_EQUAL( 1.0, bw_engine.get_lambda()); //With SI lambda is 1
}

//Test get/set lambda for breit_wheeler_engine (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_single_1 )
{
    stl_rng_wrapper rng{390109317};

    float lambda = 800.0f*flt_si_nanometer;
    breit_wheeler_engine<float, stl_rng_wrapper> bw_engine(std::move(rng));

    bw_engine.set_lambda(lambda);
    BOOST_CHECK_EQUAL( 1.0, bw_engine.get_lambda()); //With SI lambda is 1
}
