//####### Test module for picsar_span ####################################

//Define Module name
 #define BOOST_TEST_MODULE "containers/picsar_span"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <array>
#include <vector>

#include "picsar_span.hpp"

using namespace picsar::multi_physics::containers;

// ------------- Tests --------------

//***Test empty constructor

BOOST_AUTO_TEST_CASE( picsar_span_empty_constructor )
{
    auto span = picsar_span<int>();
    BOOST_CHECK_EQUAL(span.size(), 0);
    BOOST_CHECK_EQUAL(span.data(), nullptr);
}

//*******************************

//***Test constructor with raw pointers

BOOST_AUTO_TEST_CASE( picsar_span_raw_pointers_constructor )
{
    auto arr = picsar_array<double,3>{1.0,2.0,3.0};
    const auto span = picsar_span<const double>{arr.size(), arr.data()};

    BOOST_CHECK_EQUAL(span.size(), arr.size());
    for(size_t i = 0; i <  arr.size(); ++i)
        BOOST_CHECK_EQUAL(span[i], arr[i]);
}

//*******************************

//***Test range based loops

BOOST_AUTO_TEST_CASE( picsar_span_range_based_loops )
{
    auto arr = picsar_array<int,3>();
    int i = 0;
    for(auto& el : arr) el = ++i;
    const auto span = picsar_span<const int>{arr.size(), arr.data()};
    int sum = 0;
    for(const auto& el : span) sum += el;
    for(auto el : span) sum += el;
    BOOST_CHECK_EQUAL(sum, 12);
}

//*******************************

//***Test copy

BOOST_AUTO_TEST_CASE( picsar_span_copy )
{
    auto arr = picsar_array<double,3>{1.0,2.0,3.0};
    const auto span = picsar_span<double>{arr.size(), arr.data()};
    const auto cspan = span;

    BOOST_CHECK_EQUAL(cspan.size(), arr.size());
    BOOST_CHECK_EQUAL(cspan.data(), arr.data());
    for(size_t i = 0; i < arr.size(); ++i)
        BOOST_CHECK_EQUAL(cspan[i], arr[i]);
}

//*******************************
