//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/messages"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>

#include<sstream>

#include "msg.hpp"

using namespace picsar::multi_physics::utils;

// ------------- Tests --------------

//Test msg to output stream
BOOST_AUTO_TEST_CASE( msg_hello_world )
{
    std::stringbuf buf;
    std::ostream stream(&buf);
    msg("Hello world!", &stream);
    BOOST_CHECK_EQUAL( buf.str(), "Hello world!");
}

