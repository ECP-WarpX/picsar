//####### Test module for msg ####################################

//Define Module name
 #define BOOST_TEST_MODULE "chi functions"

//Will automatically define a main for this test
#define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#define PXRMP_USE_SI_UNITS
#include "chi_functions.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-5;

//Test chi functions with double precision
BOOST_AUTO_TEST_CASE( t1 )
{

    BOOST_CHECK_SMALL(1.0 , double_tolerance);
}
