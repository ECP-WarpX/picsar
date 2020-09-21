//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "math/vec_functions"

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "vec_functions.hpp"

using namespace picsar::multi_physics::math;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

// ------------- Tests --------------

// ***Test norm_square

template<typename RealType>
void test_norm_square()
{
    const auto c0 = static_cast<RealType>(1.0);
    const auto c1 = static_cast<RealType>(2.0);
    const auto c2 = static_cast<RealType>(3.0);
    const auto expected = c0*c0 + c1*c1 + c2*c2;

    const auto vec = vec3<RealType>{c0, c1, c2};

    BOOST_CHECK_SMALL((norm_square(vec)-expected)/expected, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_norm_square )
{
    test_norm_square<double>();
    test_norm_square<float>();
}

// *******************************

// ***Test norm

template<typename RealType>
void test_norm()
{
    const auto c0 = static_cast<RealType>(1.0);
    const auto c1 = static_cast<RealType>(2.0);
    const auto c2 = static_cast<RealType>(3.0);
    const auto expected = static_cast<RealType>(sqrt(c0*c0 + c1*c1 + c2*c2));

    const auto vec = vec3<RealType>{c0, c1, c2};

    BOOST_CHECK_SMALL((norm(vec)-expected)/expected, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_norm )
{
    test_norm<double>();
    test_norm<float>();
}

// *******************************

// ***Test dot product

template<typename RealType>
void test_dot()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto b0 = static_cast<RealType>(-3.0);
    const auto b1 = static_cast<RealType>(-5.0);
    const auto b2 = static_cast<RealType>(-8.0);
    const auto expected = a0*b0 + a1*b1 + a2*b2;

    const auto veca = vec3<RealType>{a0, a1, a2};
    const auto vecb = vec3<RealType>{b0, b1, b2};

    BOOST_CHECK_SMALL((dot(veca,vecb)-expected)/expected, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_dot )
{
    test_dot<double>();
    test_dot<float>();
}

// *******************************

// ***Test cross product

template<typename RealType>
void test_cross()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto b0 = static_cast<RealType>(-3.0);
    const auto b1 = static_cast<RealType>(-5.0);
    const auto b2 = static_cast<RealType>(-8.0);
    const auto exp0 = a1*b2 - a2*b1;
    const auto exp1 = a2*b0 - a0*b2;
    const auto exp2 = a0*b1 - a1*b0;

    const auto veca = vec3<RealType>{a0, a1, a2};
    const auto vecb = vec3<RealType>{b0, b1, b2};
    const auto veca_x_vecb = cross(veca, vecb);
    const auto expected = std::array<RealType,3>{exp0, exp1, exp2};

    for(int i = 0; i < 3; i++){
        if(expected[i] != 0.0)
            BOOST_CHECK_SMALL((veca_x_vecb[i]-expected[i])/expected[i], tolerance<RealType>());
        else
            BOOST_CHECK_SMALL(veca_x_vecb[i]-expected[i], tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_cross )
{
    test_cross<double>();
    test_cross<float>();
}

// *******************************

// ***Test scalar product

template<typename RealType>
void test_vec_times_scalar()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(-2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto c = static_cast<RealType>(-0.5);

    const auto exp0 = a0*c;
    const auto exp1 = a1*c;
    const auto exp2 = a2*c;

    const auto veca = vec3<RealType>{a0, a1, a2};

    const auto v_times_c = veca*c;
    const auto exptected = std::array<RealType,3>{exp0, exp1, exp2};

    for(int i = 0; i < 3; i++){
        if(exptected[i] != 0.0)
            BOOST_CHECK_SMALL((v_times_c[i]-exptected[i])/exptected[i], tolerance<RealType>());
        else
            BOOST_CHECK_SMALL(v_times_c[i]-exptected[i], tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_vec_times_scalar )
{
    test_vec_times_scalar<double>();
    test_vec_times_scalar<float>();
}

template<typename RealType>
void test_scalar_times_vec()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(-2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto c = static_cast<RealType>(-0.5);

    const auto exp0 = a0*c;
    const auto exp1 = a1*c;
    const auto exp2 = a2*c;

    const auto veca = vec3<RealType>{a0, a1, a2};

    const auto c_times_v = c*veca;
    const auto exptected = std::array<RealType,3>{exp0, exp1, exp2};

    for(int i = 0; i < 3; i++){
        if(exptected[i] != 0.0)
            BOOST_CHECK_SMALL((c_times_v[i]-exptected[i])/exptected[i], tolerance<RealType>());
        else
            BOOST_CHECK_SMALL(c_times_v[i]-exptected[i], tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_scalar_times_vec )
{
    test_scalar_times_vec<double>();
    test_scalar_times_vec<float>();
}

// *******************************

// ***Test division by scalar

template<typename RealType>
void test_vec_div_scalar()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(-2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto c = static_cast<RealType>(-0.5);

    const auto exp0 = a0/c;
    const auto exp1 = a1/c;
    const auto exp2 = a2/c;

    const auto veca = vec3<RealType>{a0, a1, a2};

    const auto v_div_c = veca/c;
    const auto exptected = std::array<RealType,3>{exp0, exp1, exp2};

    for(int i = 0; i < 3; i++){
        if(exptected[i] != 0.0)
            BOOST_CHECK_SMALL((v_div_c[i]-exptected[i])/exptected[i], tolerance<RealType>());
        else
            BOOST_CHECK_SMALL(v_div_c[i]-exptected[i], tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_vec_div_scalar )
{
    test_vec_div_scalar<double>();
    test_vec_div_scalar<float>();
}

// *******************************

// ***Test vector addition

template<typename RealType>
void test_add()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto b0 = static_cast<RealType>(-3.0);
    const auto b1 = static_cast<RealType>(-5.0);
    const auto b2 = static_cast<RealType>(-8.0);
    const auto exp0 = a0 + b0;
    const auto exp1 = a1 + b1;
    const auto exp2 = a2 + b2;

    const auto veca = vec3<RealType>{a0, a1, a2};
    const auto vecb = vec3<RealType>{b0, b1, b2};
    const auto veca_p_vecb = veca + vecb;
    const auto expected = std::array<RealType,3>{exp0, exp1, exp2};

    for(int i = 0; i < 3; i++){
        if(expected[i] != 0.0)
            BOOST_CHECK_SMALL((veca_p_vecb[i]-expected[i])/expected[i], tolerance<RealType>());
        else
            BOOST_CHECK_SMALL(veca_p_vecb[i]-expected[i], tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_add )
{
    test_add<double>();
    test_add<float>();
}

// *******************************

// ***Test vector subraction

template<typename RealType>
void test_subtract()
{
    const auto a0 = static_cast<RealType>(1.0);
    const auto a1 = static_cast<RealType>(2.0);
    const auto a2 = static_cast<RealType>(3.0);
    const auto b0 = static_cast<RealType>(-3.0);
    const auto b1 = static_cast<RealType>(-5.0);
    const auto b2 = static_cast<RealType>(-8.0);
    const auto exp0 = a0 - b0;
    const auto exp1 = a1 - b1;
    const auto exp2 = a2 - b2;

    const auto veca = vec3<RealType>{a0, a1, a2};
    const auto vecb = vec3<RealType>{b0, b1, b2};
    const auto veca_m_vecb = veca - vecb;
    const auto expected = std::array<RealType,3>{exp0, exp1, exp2};

    for(int i = 0; i < 3; i++){
        if(expected[i] != 0.0)
            BOOST_CHECK_SMALL((veca_m_vecb[i]-expected[i])/expected[i], tolerance<RealType>());
        else
            BOOST_CHECK_SMALL(veca_m_vecb[i]-expected[i], tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_vec_functions_subtract )
{
    test_subtract<double>();
    test_subtract<float>();
}

// *******************************
