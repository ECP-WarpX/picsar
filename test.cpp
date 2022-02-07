#define BOOST_TEST_MODULE Test
#define BOOST_ALL_NO_LIB
#define BOOST_UNIT_TEST_FRAMEWORK_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test )
{
    boost::test_tools::output_test_stream output;
    const int a = 1;
    output << a;
    BOOST_CHECK(output.is_equal("1"));
}
