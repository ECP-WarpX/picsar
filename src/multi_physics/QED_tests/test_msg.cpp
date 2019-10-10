//####### Test module for msg ####################################

//Define Module name
 #define BOOST_TEST_MODULE "messages"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

#include<sstream>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS

#include "msg.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Test msg to output stream
BOOST_AUTO_TEST_CASE( msg_hello_world )
{
    std::stringbuf buf;
    std::ostream stream(&buf);
    msg("Hello world!", &stream);
    BOOST_CHECK_EQUAL( buf.str(), "Hello world!");
}
