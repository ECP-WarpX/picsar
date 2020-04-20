//####### Test module for picsar_algo ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/picsar_algo"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <array>
#include <algorithm>

//Force the use of PICSAR implementation of upper_bound for debug purposes
#define PXRMP_FORCE_PICSAR_UPPER_BOUND
#include "picsar_algo.hpp"

using namespace picsar::multi_physics::utils;

// ------------- Tests --------------

// ***Test upper_bound

BOOST_AUTO_TEST_CASE( picsar_upper_bound_1 )
{
    const auto arr = std::array<double,5>{0.0, 1.0, 2.0, 3.0, 4.0};

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), -1.0),
        std::upper_bound(arr.begin(), arr.end(), -1.0));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 0.0),
        std::upper_bound(arr.begin(), arr.end(), 0.0));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 0.1),
        std::upper_bound(arr.begin(), arr.end(), 0.1));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 1.0),
        std::upper_bound(arr.begin(), arr.end(), 1.0));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 1.1),
        std::upper_bound(arr.begin(), arr.end(), 1.1));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 3.9),
        std::upper_bound(arr.begin(), arr.end(), 3.9));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 4.0),
        std::upper_bound(arr.begin(), arr.end(), 4.0));

    BOOST_CHECK_EQUAL(
        picsar_upper_bound(arr.begin(), arr.end(), 5.0),
        std::upper_bound(arr.begin(), arr.end(), 5.0));
}

// *******************************
