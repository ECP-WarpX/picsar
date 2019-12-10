//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "containers/picsar_array"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#define PXRMP_FORCE_PICSAR_ARRAY
#include "picsar_array.hpp"

using namespace picsar::multi_physics::containers;


// ------------- Tests --------------

//Test empty constructor
BOOST_AUTO_TEST_CASE( picsar_array_empty_constructor )
{
    auto arr = picsar_array<int,3>();
    arr[0] = 3;
    arr[1] = 2;
    arr[2] = 1;
    BOOST_CHECK_EQUAL(arr[0], 3);
    BOOST_CHECK_EQUAL(arr[1], 2);
    BOOST_CHECK_EQUAL(arr[2], 1);
}

//Test initializer list
BOOST_AUTO_TEST_CASE( picsar_array_list_constructor )
{
    auto arr = picsar_array<int,3>{3,2,1};
    const auto carr = picsar_array<int,3>{3,2,1};
    BOOST_CHECK_EQUAL(arr[0], 3);
    BOOST_CHECK_EQUAL(arr[1], 2);
    BOOST_CHECK_EQUAL(arr[2], 1);
    BOOST_CHECK_EQUAL(carr[0], 3);
    BOOST_CHECK_EQUAL(carr[1], 2);
    BOOST_CHECK_EQUAL(carr[2], 1);
}

//Test iterator
BOOST_AUTO_TEST_CASE( picsar_array_range_based_loops )
{
    auto arr = picsar_array<int,3>();
    int i = 0;
    for(auto& el : arr) el = ++i;
    int sum = 0;
    for(const auto& el : arr) sum += el;
    for(auto el : arr) sum += el;
    BOOST_CHECK_EQUAL(sum, 12);
}
