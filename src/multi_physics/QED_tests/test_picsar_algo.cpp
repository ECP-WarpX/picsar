//####### Test module for picsar_algo ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/picsar_algo"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

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

    for (const auto nn : std::vector<double>{-1.,0.,0.1,1,1.1,3.9,4.0,5.0}){
        BOOST_CHECK_EQUAL(
            picsar_upper_bound(arr.begin(), arr.end(), nn),
            std::upper_bound(arr.begin(), arr.end(), nn));
    }
}

// ***Test upper_bound_functor

BOOST_AUTO_TEST_CASE( picsar_upper_bound_functor_1 )
{
    const auto arr = std::array<double,5>{0.0, 1.0, 2.0, 3.0, 4.0};

    for (const auto nn : std::vector<double>{-1.,0.,0.1,1,1.1,3.9,4.0,5.0}){
        BOOST_CHECK_EQUAL(
            picsar_upper_bound_functor(0, arr.size(), nn, [&](size_t i){
                return arr[i];}),
            std::distance(
                arr.begin(),
                std::upper_bound(arr.begin(), arr.end(), nn)));

    }
}

// *******************************
