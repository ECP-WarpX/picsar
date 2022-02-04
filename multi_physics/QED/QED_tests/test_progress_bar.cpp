//####### Test module for progress bar ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/progress_bar"

//Include Boost unit tests library & library for out stream tests
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

// ------------- Tests --------------

// ***Test progress bar

BOOST_AUTO_TEST_CASE( picsar_progress_bar_1 )
{
    boost::test_tools::output_test_stream output;
    const int a = 1;
    output << a;
    BOOST_CHECK(output.is_equal("1"));
}

// *******************************
