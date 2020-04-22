//####### Test module for chi functions #############################

//Define Module name
 #define BOOST_TEST_MODULE "phys/chi_functions"

//Will automatically define a main for this test
#define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "chi_functions.hpp"

#include "unit_conversion.hpp"
#include "vec_functions.hpp"

using namespace picsar::multi_physics::phys;
using namespace picsar::multi_physics::math;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-6;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}


const double me_c = electron_mass<double> * light_speed<double>;
const double lambda = 800.0e-9;
const double omega = 2.0*pi<double>*light_speed<double>/lambda;
const double eref = 2.0*pi<double>*electron_mass<double>*
    light_speed<double>*light_speed<double>/(lambda*elementary_charge<double>);
const double bref = eref/light_speed<double>;

// ------------- Tests --------------

//***Test chi calculation for photons

template<typename RealType, unit_system UnitSystem>
void test_chi_photons(
    vec3<double> t_p, vec3<double> t_em_e, vec3<double> t_em_b,
    double t_exp_res, double t_ref_quant = 1.0)
{
    const auto p =  vec3<RealType>{
        static_cast<RealType>(t_p[0]),
        static_cast<RealType>(t_p[1]),
        static_cast<RealType>(t_p[2])}*
        conv<quantity::momentum, unit_system::SI, UnitSystem, RealType>::fact();

    const auto em_e =  vec3<RealType>{
        static_cast<RealType>(t_em_e[0]),
        static_cast<RealType>(t_em_e[1]),
        static_cast<RealType>(t_em_e[2])}*
        conv<quantity::E, unit_system::SI, UnitSystem, RealType>::fact(1.0, t_ref_quant);

    const auto em_b =  vec3<RealType>{
        static_cast<RealType>(t_em_b[0]),
        static_cast<RealType>(t_em_b[1]),
        static_cast<RealType>(t_em_b[2])}*
        conv<quantity::B, unit_system::SI, UnitSystem, RealType>::fact(1.0, t_ref_quant);

    const auto res = chi_photon<RealType, UnitSystem>(p, em_e, em_b, t_ref_quant);

    const auto exp_res = static_cast<RealType>(t_exp_res);

    if(exp_res != static_cast<RealType>(0.0))
        BOOST_CHECK_SMALL((res-exp_res)/exp_res, tolerance<RealType>());
    else
        BOOST_CHECK_SMALL(res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( chi_photons_1 )
{
    const double px = 83.759*me_c;
    const double py = 139.311*me_c;
    const double pz = -230.553*me_c;
    const double ex = -166.145*eref;
    const double ey = -78.231*eref;
    const double ez = -278.856*eref;
    const double bx = -279.174*bref;
    const double by = -158.849*bref;
    const double bz = -93.826*bref;

    const double chi_exp = 0.3471118445206898;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_photons<double, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_photons<float, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_photons_2 )
{
    const double px = 9.2627*me_c;
    const double py = -25.4575*me_c;
    const double pz = -10.2246*me_c;
    const double ex = 2.9271*eref;
    const double ey = 10.4293*eref;
    const double ez = 3.6103*eref;
    const double bx = 1.7439*bref;
    const double by = 1.9778*bref;
    const double bz = 17.8799*bref;

    const double chi_exp = 0.0009041474058668584;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_photons<double, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_photons<float, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_photons_3 )
{
    const double px = -2314.45*me_c;
    const double py = -2356.30*me_c;
    const double pz = 546.28*me_c;
    const double ex = 1230.11*eref;
    const double ey = 1638.02*eref;
    const double ez = -2911.04*eref;
    const double bx = -2203.66*bref;
    const double by = 1243.79*bref;
    const double bz = -2830.99*bref;

    const double chi_exp = 57.22043983047014;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_photons<double, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_photons<float, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_photons_4 )
{
    const double px = 0*me_c;
    const double py = 0*me_c;
    const double pz = 0*me_c;
    const double ex = 1230.11*eref;
    const double ey = 1638.02*eref;
    const double ez = -2911.04*eref;
    const double bx = -2203.66*bref;
    const double by = 1243.79*bref;
    const double bz = -2830.99*bref;

    const double chi_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_photons<double, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_photons<float, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_photons_5 )
{
    const double px = -2314.45*me_c;
    const double py = -2356.30*me_c;
    const double pz = 546.28*me_c;
    const double ex = 0.0*eref;
    const double ey = 0.0*eref;
    const double ez = 0.0*eref;
    const double bx = 0.0*bref;
    const double by = 0.0*bref;
    const double bz = 0.0*bref;

    const double chi_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_photons<double, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_photons<float, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_photons_6 )
{
    const double px = -2314.45*me_c;
    const double py = 0.0*me_c;
    const double pz = 0.0*me_c;
    const double ex = 1230.11*eref;
    const double ey = 0.0*eref;
    const double ez = 0.0*eref;
    const double bx = 0.0*bref;
    const double by = 0.0*bref;
    const double bz = 0.0*bref;

    const double chi_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_photons<double, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_photons<float, unit_system::SI>(p, em_e, em_b, chi_exp);
    test_chi_photons<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_photons<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_photons<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

//*******************************

//***Test chi calculation for electrons and positrons

template<typename RealType, unit_system UnitSystem>
void test_chi_ele_pos(
    vec3<double> t_p, vec3<double> t_em_e, vec3<double> t_em_b,
    double t_exp_res, double t_ref_quant = 1.0)
{
    const auto p =  vec3<RealType>{
        static_cast<RealType>(t_p[0]),
        static_cast<RealType>(t_p[1]),
        static_cast<RealType>(t_p[2])}*
        conv<quantity::momentum, unit_system::SI, UnitSystem, RealType>::fact();

    const auto em_e =  vec3<RealType>{
        static_cast<RealType>(t_em_e[0]),
        static_cast<RealType>(t_em_e[1]),
        static_cast<RealType>(t_em_e[2])}*
        conv<quantity::E, unit_system::SI, UnitSystem, RealType>::fact(1.0, t_ref_quant);

    const auto em_b =  vec3<RealType>{
        static_cast<RealType>(t_em_b[0]),
        static_cast<RealType>(t_em_b[1]),
        static_cast<RealType>(t_em_b[2])}*
        conv<quantity::B, unit_system::SI, UnitSystem, RealType>::fact(1.0, t_ref_quant);

    const auto res = chi_ele_pos<RealType, UnitSystem>(p, em_e, em_b, t_ref_quant);

    const auto exp_res = static_cast<RealType>(t_exp_res);

    if(exp_res != static_cast<RealType>(0.0))
        BOOST_CHECK_SMALL((res-exp_res)/exp_res, tolerance<RealType>());
    else
        BOOST_CHECK_SMALL(res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( chi_ele_pos_1 )
{
    const double px = 24.3752*me_c;
    const double py = -11.5710*me_c;
    const double pz = -10.0841*me_c;
    const double ex = 57.185*eref;
    const double ey = -16.6555*eref;
    const double ez = 22.4340*eref;
    const double bx = 6.6911*bref;
    const double by = -23.8724*bref;
    const double bz = 13.9934*bref;

    const double chi_exp = 0.002167166273468506;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_ele_pos<double, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);

}

BOOST_AUTO_TEST_CASE( chi_ele_pos_2 )
{
    const double px = 4.015*me_c;
    const double py = 197.287*me_c;
    const double pz = 141.705*me_c;
    const double ex = 30.287*eref;
    const double ey = 115.740*eref;
    const double ez = 120.891*eref;
    const double bx = -190.161*bref;
    const double by = -129.115*bref;
    const double bz = -57.002*bref;

    const double  chi_exp = 0.16631811297207244;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_ele_pos<double, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_ele_pos_3 )
{
    const double px = -2534.83*me_c;
    const double py = 1011.54*me_c;
    const double pz = -793.04*me_c;
    const double ex = 741.67*eref;
    const double ey = -2359.97*eref;
    const double ez = 1463.50*eref;
    const double bx = 1477.19*bref;
    const double by = -1448.33*bref;
    const double bz = 1953.68*bref;

    const double chi_exp = 16.011457274095676;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_ele_pos<double, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_ele_pos_4 )
{
    const double px = 0*me_c;
    const double py = 0*me_c;
    const double pz = 0*me_c;
    const double ex = 741.67*eref;
    const double ey = -2359.97*eref;
    const double ez = 1463.50*eref;
    const double bx = 1477.19*bref;
    const double by = -1448.33*bref;
    const double bz = 1953.68*bref;

    const double chi_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_ele_pos<double, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_ele_pos_5 )
{
    const double px = -2534.83*me_c;
    const double py = 1011.54*me_c;
    const double pz = -793.04*me_c;
    const double ex = 0*eref;
    const double ey = 0*eref;
    const double ez = 0*eref;
    const double bx = 0*bref;
    const double by = 0*bref;
    const double bz = 0*bref;

    const double chi_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_ele_pos<double, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( chi_ele_pos_6 )
{
    const double px = -2534.83*me_c;
    const double py = 0*me_c;
    const double pz = 0*me_c;
    const double ex = 0*eref;
    const double ey = 0*eref;
    const double ez = 0*eref;
    const double bx = 1477.19*bref;
    const double by = 0*bref;
    const double bz = 0*bref;

    const double chi_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};
    const auto em_e = vec3<double>{ex,ey,ez};
    const auto em_b = vec3<double>{bx,by,bz};

    test_chi_ele_pos<double, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<double, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<double, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<double, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::SI>(p, em_e, em_b, chi_exp, 1.0);
    test_chi_ele_pos<float, unit_system::norm_lambda>(p, em_e, em_b, chi_exp, lambda);
    test_chi_ele_pos<float, unit_system::norm_omega>(p, em_e, em_b, chi_exp, omega);
    test_chi_ele_pos<float, unit_system::heaviside_lorentz>(p, em_e, em_b, chi_exp, 1.0);
}

//*******************************
