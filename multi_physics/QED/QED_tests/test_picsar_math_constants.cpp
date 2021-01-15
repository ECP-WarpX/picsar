//####### Test module for mathematical constants ##################

//Define Module name
 #define BOOST_TEST_MODULE "math/constants"

#include <cmath>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <picsar_qed/math/math_constants.h>

using namespace picsar::multi_physics::math;

// ------------- Tests --------------

// ***Test math constants

template<typename RealType>
void test_case_const_math()
{
    const auto exp_pi =
        static_cast<RealType>(3.14159265358979323846264338327950288);

    const auto exp_zero = static_cast<RealType>(0.0);
    const auto exp_half = static_cast<RealType>(0.5);
    const auto exp_one = static_cast<RealType>(1.0);
    const auto exp_two = static_cast<RealType>(2.0);
    const auto exp_three = static_cast<RealType>(3.0);
    const auto exp_four = static_cast<RealType>(4.0);
    const auto exp_one_third = static_cast<RealType>(1.0/3.0);
    const auto exp_two_thirds = static_cast<RealType>(2.0/3.0);
    const auto exp_five_thirds = static_cast<RealType>(5.0/3.0);

    BOOST_CHECK_EQUAL(pi<RealType>, exp_pi);
    BOOST_CHECK_EQUAL(zero<RealType>, exp_zero);
    BOOST_CHECK_EQUAL(half<RealType>, exp_half);
    BOOST_CHECK_EQUAL(one<RealType>, exp_one);
    BOOST_CHECK_EQUAL(two<RealType>, exp_two);
    BOOST_CHECK_EQUAL(three<RealType>, exp_three);
    BOOST_CHECK_EQUAL(four<RealType>, exp_four);
    BOOST_CHECK_EQUAL(one_third<RealType>, exp_one_third);
    BOOST_CHECK_EQUAL(two_thirds<RealType>, exp_two_thirds);
    BOOST_CHECK_EQUAL(five_thirds<RealType>, exp_five_thirds);
}

BOOST_AUTO_TEST_CASE( picsar_const_math )
{
    test_case_const_math<double>();
    test_case_const_math<float>();
}

// *******************************
