//####### Test module for utilites ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utilites"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 //Include Boost unit tests library & library for floating point comparison
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS

#include "utilities.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-10;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//Templated tolerance
template <typename T>
T tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//Test generate_lin_spaced_vec generic
template<typename _WHATEVER>
void test_generate_lin_spaced_vec(_WHATEVER min, _WHATEVER max, size_t size)
{
    std::vector<_WHATEVER> vv =
        generate_lin_spaced_vec<_WHATEVER>(min, max, size);

    BOOST_CHECK_EQUAL(vv.size(), size);
    BOOST_CHECK_EQUAL(vv.front(), min);
    BOOST_CHECK_EQUAL(vv.back(), max);

    _WHATEVER expdiff =(max-min)/(size-1) ;

    for(size_t i = 1; i < vv.size(); ++i){
        _WHATEVER diff = vv[i]-vv[i-1];
        BOOST_CHECK_SMALL((diff-expdiff)/expdiff, tolerance<_WHATEVER>());
    }
}


//Test generate_lin_spaced_vec generic with double precision
BOOST_AUTO_TEST_CASE( test_generate_lin_spaced_vec_double_1 )
{
    test_generate_lin_spaced_vec<double>(-73, 112, 500);
}

//Test generate_lin_spaced_vec generic with single precision
BOOST_AUTO_TEST_CASE( test_generate_lin_spaced_vec_single_1 )
{
    test_generate_lin_spaced_vec<float>(-73, 112, 500);
}


//Test generate_log_spaced_vec generic
template<typename _WHATEVER>
void test_generate_log_spaced_vec(_WHATEVER min, _WHATEVER max, size_t size)
{
    std::vector<_WHATEVER> vv =
        generate_log_spaced_vec<_WHATEVER>(min, max, size);

    BOOST_CHECK_EQUAL(vv.size(), size);
    BOOST_CHECK_EQUAL(vv.front(), min);
    BOOST_CHECK_EQUAL(vv.back(), max);

    _WHATEVER mul = pow(max/min, static_cast<_WHATEVER>(1.0/(size-1)));

    for(size_t i = 1; i < vv.size(); ++i){
        _WHATEVER ratio = vv[i]/vv[i-1];
        BOOST_CHECK_SMALL((ratio-mul)/mul, tolerance<_WHATEVER>());
    }
}


//Test generate_log_spaced_vec generic with double precision
BOOST_AUTO_TEST_CASE( test_generate_log_spaced_vec_double_1 )
{
    test_generate_log_spaced_vec<double>(0.01, 100.0, 200);
}

//Test generate_log_spaced_vec generic with single precision
BOOST_AUTO_TEST_CASE( test_generate_log_spaced_vec_single_1 )
{
    test_generate_log_spaced_vec<float>(0.01, 100.0, 200);
}
