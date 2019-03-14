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

// ------------- Tests --------------

//Test get/set lambda for breit_wheeler_engine (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_double_1 )
{
  //Cumbersome, but C++11 does not have make_unique :-(
  auto ptr_wrp = std::unique_ptr<rng_wrapper<double>>
  (new stl_rng_wrapper<double>(390109317));

  double lambda = 800.0*si_nanometer;
  breit_wheeler_engine<double> bw_engine(move(ptr_wrp));
  bw_engine.set_lambda(lambda);
  BOOST_CHECK_EQUAL( lambda, bw_engine.get_lambda());
}

//Test get/set lambda for breit_wheeler_engine (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_single_1 )
{
  //Cumbersome, but C++11 does not have make_unique :-(
  auto ptr_wrp = std::unique_ptr<rng_wrapper<float>>
  (new stl_rng_wrapper<float>(390109317));

  float lambda = 800.0f*flt_si_nanometer;
  breit_wheeler_engine<float> bw_engine(move(ptr_wrp));
  bw_engine.set_lambda(lambda);
  BOOST_CHECK_EQUAL( lambda, bw_engine.get_lambda());
}
