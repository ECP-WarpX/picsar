//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "vector functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS

#include "vec_functions.hpp"

using namespace picsar::multi_physics;

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

// ------------- Tests --------------

//Test norm2 generic
template<typename T>
void vec_functions_norm2()
{
    vec3<T> vv{static_cast<T>(1.0),static_cast<T>(-2.0),static_cast<T>(3.0)};
    T exp = static_cast<T>(1.0 + 4.0 + 9.0);
    BOOST_CHECK_SMALL((norm2(vv)-exp)/exp, tolerance<T>());
}

//Test norm2 in double precision
BOOST_AUTO_TEST_CASE( vec_functions_norm2_double_1 )
{
    vec_functions_norm2<double>();
}

//Test norm2 in single precision
BOOST_AUTO_TEST_CASE( vec_functions_norm2_single_1 )
{
    vec_functions_norm2<float>();
}

//Test norm generic
template<typename T>
void vec_functions_norm()
{
    vec3<T> vv{static_cast<T>(1.0),static_cast<T>(-2.0),static_cast<T>(3.0)};
    T exp =  static_cast<T>(sqrt(1.0 + 4.0 + 9.0));
    BOOST_CHECK_SMALL((norm(vv)-exp)/exp,  tolerance<T>());
}

//Test norm in double precision
BOOST_AUTO_TEST_CASE( vec_functions_norm_double_1 )
{
    vec_functions_norm<double>();
}

//Test norm in single precision
BOOST_AUTO_TEST_CASE( vec_functions_norm_single_1 )
{
    vec_functions_norm<float>();
}

//Test dot generic
template<typename T>
void vec_functions_dot()
{
    vec3<T> vv1{static_cast<T>(1.0),static_cast<T>(-2.0),static_cast<T>(3.0)};
    vec3<T> vv2{static_cast<T>(-1.0),static_cast<T>(0.0),static_cast<T>(5.0)};
    T exp = static_cast<T>(14.0);
    BOOST_CHECK_SMALL((dot(vv1, vv2)-exp)/exp, tolerance<T>());
}

//Test dot in double precision
BOOST_AUTO_TEST_CASE( vec_functions_dot_double_1 )
{
    vec_functions_dot<double>();
}

//Test dot in single precision
BOOST_AUTO_TEST_CASE( vec_functions_dot_single_1 )
{
    vec_functions_dot<float>();
}

//Test cross generic
template<typename T>
void vec_functions_cross()
{
    vec3<T> vv1{static_cast<T>(1./3.),static_cast<T>(-1./4.),static_cast<T>(1./5.)};
    vec3<T> vv2{static_cast<T>(-2./3.),static_cast<T>(-3./4.),static_cast<T>(4./5.)};
    vec3<T> exp{static_cast<T>(-1./20.),static_cast<T>(-2./5.),static_cast<T>(-5./12.)};

    vec3<T> res = cross(vv1,vv2);
    BOOST_CHECK_SMALL((res[0]-exp[0])/exp[0], tolerance<T>());
    BOOST_CHECK_SMALL((res[1]-exp[1])/exp[1], tolerance<T>());
    BOOST_CHECK_SMALL((res[2]-exp[2])/exp[2], tolerance<T>());
}

//Test cross in double precision
BOOST_AUTO_TEST_CASE( vec_functions_cross_double_1 )
{
    vec_functions_cross<double>();
}

//Test cross in single precision
BOOST_AUTO_TEST_CASE( vec_functions_cross_single_1 )
{
    vec_functions_cross<float>();
}

//Test vector times scalar generic
template<typename T>
void vec_functions_vsprod()
{
    vec3<T> vv{static_cast<T>(1.0),static_cast<T>(2.0),static_cast<T>(-3.0)};
    T s = static_cast<T>(-2.0);
    vec3<T> exp{static_cast<T>(-2.0),static_cast<T>(-4.0),static_cast<T>(6.0)};

    vec3<T> r1 = s * vv;
    vec3<T> r2 = vv * s;

    BOOST_CHECK_SMALL((r1[0]-exp[0])/exp[0], tolerance<T>());
    BOOST_CHECK_SMALL((r1[1]-exp[1])/exp[1], tolerance<T>());
    BOOST_CHECK_SMALL((r1[2]-exp[2])/exp[2], tolerance<T>());
    BOOST_CHECK_SMALL((r2[0]-exp[0])/exp[0], tolerance<T>());
    BOOST_CHECK_SMALL((r2[1]-exp[1])/exp[1], tolerance<T>());
    BOOST_CHECK_SMALL((r2[2]-exp[2])/exp[2], tolerance<T>());
}

//Test vector times scalar in double precision
BOOST_AUTO_TEST_CASE( vec_functions_vsprod_double_1 )
{
    vec_functions_vsprod<double>();
}

//Test vector times scalar in single precision
BOOST_AUTO_TEST_CASE( vec_functions_vsprod_float_1 )
{
    vec_functions_vsprod<float>();
}

//Test vector divided by scalar generic
template<typename T>
void vec_functions_vsdiv()
{
    vec3<T> vv{static_cast<T>(2.0),static_cast<T>(4.0),static_cast<T>(-6.0)};
    T s = static_cast<T>(-2.0);
    vec3<T> exp{static_cast<T>(-1.0),static_cast<T>(-2.0),static_cast<T>(3.0)};

    vec3<T> r = vv / s;

    BOOST_CHECK_SMALL((r[0]-exp[0])/exp[0], tolerance<T>());
    BOOST_CHECK_SMALL((r[1]-exp[1])/exp[1], tolerance<T>());
    BOOST_CHECK_SMALL((r[2]-exp[2])/exp[2], tolerance<T>());
}

//Test vector divided by scalar in double precision
BOOST_AUTO_TEST_CASE( vec_functions_vsdiv_double_1 )
{
    vec_functions_vsdiv<double>();
}

//Test vector divided by scalar in single precision
BOOST_AUTO_TEST_CASE( vec_functions_vsdiv_float_1 )
{
    vec_functions_vsdiv<float>();
}

//Test vector add generic
template<typename T>
void vec_functions_vadd()
{
    vec3<T> vv1{static_cast<T>(1.0),static_cast<T>(2.0),static_cast<T>(3.0)};
    vec3<T> vv2{static_cast<T>(1.0),static_cast<T>(-1.0),static_cast<T>(-2.0)};

    vec3<T> exp{static_cast<T>(2.0),static_cast<T>(1.0),static_cast<T>(1.0)};

    vec3<T> r = vv1 + vv2;

    BOOST_CHECK_SMALL((exp[0]-r[0])/exp[0], tolerance<T>());
    BOOST_CHECK_SMALL((exp[1]-r[1])/exp[1], tolerance<T>());
    BOOST_CHECK_SMALL((exp[2]-r[2])/exp[2], tolerance<T>());
}

//Test vector add in double precision
BOOST_AUTO_TEST_CASE( vec_functions_vadd_double_1 )
{
    vec_functions_vadd<double>();
}

//Test vector add in single precision
BOOST_AUTO_TEST_CASE( vec_functions_vadd_single_1 )
{
    vec_functions_vadd<float>();
}

//Test vector diff generic
template<typename T>
void vec_functions_vdiff()
{
    vec3<T> vv1{static_cast<T>(1.0),static_cast<T>(2.0),static_cast<T>(3.0)};
    vec3<T> vv2{static_cast<T>(2.0),static_cast<T>(-1.0),static_cast<T>(-2.0)};

    vec3<T> exp{static_cast<T>(-1.0),static_cast<T>(3.0),static_cast<T>(5.0)};

    vec3<T> r = vv1 - vv2;

    BOOST_CHECK_SMALL((r[0]-exp[0])/exp[0], tolerance<T>());
    BOOST_CHECK_SMALL((r[1]-exp[1])/exp[1], tolerance<T>());
    BOOST_CHECK_SMALL((r[2]-exp[2])/exp[2], tolerance<T>());
}

//Test vector diff in double precision
BOOST_AUTO_TEST_CASE( vec_functions_vdiff_double_1 )
{
    vec_functions_vdiff<double>();
}

//Test vector diff in single precision
BOOST_AUTO_TEST_CASE( vec_functions_vdiff_single_1 )
{
    vec_functions_vdiff<float>();
}
