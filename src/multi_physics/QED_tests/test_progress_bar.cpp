//####### Test module for progress bar ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/progress_bar"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for out stream tests
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#include "progress_bar.hpp"

using namespace picsar::multi_physics::utils;

// ------------- Tests --------------

// ***Test progress bar

BOOST_AUTO_TEST_CASE( picsar_progress_bar_1 )
{
    boost::test_tools::output_test_stream output;
    draw_progress(0, 10, "aa", 1, false, output);
    BOOST_CHECK(output.is_equal(
        " [>                                                 ] 0%  aa\r"));
}

BOOST_AUTO_TEST_CASE( picsar_progress_bar_2 )
{
    boost::test_tools::output_test_stream output;
    draw_progress(2, 10, "bcd", 1, false, output);
    BOOST_CHECK(output.is_equal(
        " [==========>                                       ] 20%  bcd\r"));
}

BOOST_AUTO_TEST_CASE( picsar_progress_bar_3 )
{
    boost::test_tools::output_test_stream output;
    draw_progress(10, 10, "efg", 1, false, output);
    BOOST_CHECK(output.is_equal(
        " [==================================================] 100%  efg\r"));
}

BOOST_AUTO_TEST_CASE( picsar_progress_bar_4 )
{
    boost::test_tools::output_test_stream output;
    draw_progress(10, 10, "efg", 1, true, output);
    BOOST_CHECK(output.is_equal(
        " [==================================================] 100%  efg\n"));
}

BOOST_AUTO_TEST_CASE( picsar_progress_bar_5 )
{
    boost::test_tools::output_test_stream output;
    draw_progress(5, 10, "efg", 3, false, output);
    BOOST_CHECK(output.is_equal(""));
}

// *******************************
