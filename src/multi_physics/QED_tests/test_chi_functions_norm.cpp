//####### Test module for chi functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "chi functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//This module will use
#define PXRMP_USE_NORMALIZED_UNITS
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

/*
cout << "calc Chi for photons (mom=[83.759, 139.311, -230.553], EB=[-166.145, -78.231, -278.856, -279.174, -158.849, -93.826], l = 800 nm, exp. 0.347111844317) :" << endl;
cout << picsar::multi_physics::chi_photon_lambda({83.759, 139.311, -230.553},{-166.145, -78.231, -278.856, -279.174, -158.849, -93.826}, 0.8 * picsar::multi_physics::_um) << endl;
cout << "calc Chi for photons (mom=[-2314.45, -2356.30, 546.28], EB=[1230.11, 1638.02, -2911.04, -2203.66, 1243.79, -2830.99], l = 800 nm, exp. 57.2204397969) :" << endl;
cout << picsar::multi_physics::chi_photon_lambda({-2314.45, -2356.30, 546.2},{1230.11, 1638.02, -2911.04, -2203.66, 1243.79, -2830.99}, 0.8 * picsar::multi_physics::_um) << endl;
cout << "calc Chi for photons (mom=[9.2627, -25.4575, -10.2246], EB=[2.9271, 10.4293, 3.6103, 1.7439, 1.9778, 17.8799], l = 800 nm, exp. 0.000904147405336) :" << endl;
cout << picsar::multi_physics::chi_photon_lambda({9.2627, -25.4575, -10.2246},{2.9271, 10.4293, 3.6103, 1.7439, 1.9778, 17.8799}, 0.8 * picsar::multi_physics::_um) << endl;
*/
