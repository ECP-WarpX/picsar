//####### Test module for gamma functions #############################

//Define Module name
 #define BOOST_TEST_MODULE "phys/gamma_functions"

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <picsar_qed/physics/gamma_functions.hpp>

#include <picsar_qed/physics/unit_conversion.hpp>
#include <picsar_qed/math/vec_functions.hpp>

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

//***Test normalized energy calculation for photons

template<typename RealType, unit_system UnitSystem>
void test_gamma_photons(
    vec3<double> t_p, double t_exp_res, double t_ref_quant = 1.0)
{
    const auto p =  vec3<RealType>{
        static_cast<RealType>(t_p[0]),
        static_cast<RealType>(t_p[1]),
        static_cast<RealType>(t_p[2])}*
        conv<quantity::momentum, unit_system::SI, UnitSystem, RealType>::fact();

    const auto res = compute_gamma_photon<RealType, UnitSystem>(p, t_ref_quant);

    const auto exp_res = static_cast<RealType>(t_exp_res);

    if(exp_res != static_cast<RealType>(0.0))
        BOOST_CHECK_SMALL((res-exp_res)/exp_res, tolerance<RealType>());
    else
        BOOST_CHECK_SMALL(res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( gamma_photons_1 )
{
    const double px = 83.759*me_c;
    const double py = 139.311*me_c;
    const double pz = -230.553*me_c;

    const double gamma_exp = 282.09539275039566;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_photons<double, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_photons<float, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_photons_2 )
{
    const double px = 9.2627*me_c;
    const double py = -25.4575*me_c;
    const double pz = -10.2246*me_c;

    const double gamma_exp = 28.955558407670193;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_photons<double, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_photons<float, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_photons_3 )
{
    const double px = -2314.45*me_c;
    const double py = -2356.30*me_c;
    const double pz = 546.28*me_c;

    const double gamma_exp = 3347.723156251126;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_photons<double, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_photons<float, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_photons_4 )
{
    const double px = 0*me_c;
    const double py = 0*me_c;
    const double pz = 0*me_c;

    const double gamma_exp = 0.0;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_photons<double, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_photons<float, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_photons_5 )
{
    const double px = 0.001*me_c;
    const double py = -0.002*me_c;
    const double pz = 0.003*me_c;

    const double gamma_exp = 0.0037416573867739412;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_photons<double, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_photons<float, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<float, unit_system::norm_lambda>(p,gamma_exp, lambda);
    test_gamma_photons<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_photons_6 )
{
    const double px = -2314.45*me_c;
    const double py = 0.0*me_c;
    const double pz = 0.0*me_c;

    const double gamma_exp = 2314.45;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_photons<double, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_photons<float, unit_system::SI>(p, gamma_exp);
    test_gamma_photons<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_photons<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_photons<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

//*******************************

//***Test Lorentz factor calculation for electrons and positrons

template<typename RealType, unit_system UnitSystem>
void test_gamma_ele_pos(
    vec3<double> t_p, double t_exp_res, double t_ref_quant = 1.0)
{
    const auto p =  vec3<RealType>{
        static_cast<RealType>(t_p[0]),
        static_cast<RealType>(t_p[1]),
        static_cast<RealType>(t_p[2])}*
        conv<quantity::momentum, unit_system::SI, UnitSystem, RealType>::fact();

    const auto res = compute_gamma_ele_pos<RealType, UnitSystem>(p, t_ref_quant);

    const auto exp_res = static_cast<RealType>(t_exp_res);

    if(exp_res != static_cast<RealType>(0.0))
        BOOST_CHECK_SMALL((res-exp_res)/exp_res, tolerance<RealType>());
    else
        BOOST_CHECK_SMALL(res, tolerance<RealType>());
}

BOOST_AUTO_TEST_CASE( gamma_ele_pos_1 )
{
    const double px = 24.3752*me_c;
    const double py = -11.5710*me_c;
    const double pz = -10.0841*me_c;

    const double gamma_exp = 28.822343569703;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_ele_pos<double, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<float, unit_system::heaviside_lorentz>(p,gamma_exp, 1.0);

}

BOOST_AUTO_TEST_CASE( gamma_ele_pos_2 )
{
    const double px = 4.015*me_c;
    const double py = 197.287*me_c;
    const double pz = 141.705*me_c;

    const double  gamma_exp = 242.93947315946826;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_ele_pos<double, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_ele_pos_3 )
{
    const double px = -2534.83*me_c;
    const double py = 1011.54*me_c;
    const double pz = -793.04*me_c;

    const double gamma_exp = 2842.0924935863713;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_ele_pos<double, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_ele_pos_4 )
{
    const double px = 0*me_c;
    const double py = 0*me_c;
    const double pz = 0*me_c;

    const double gamma_exp = 1.0;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_ele_pos<double, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_ele_pos_5 )
{
    const double px = 0.001*me_c;
    const double py = -0.002*me_c;
    const double pz = 0.003*me_c;

    const double gamma_exp = 1.0000069999755001;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_ele_pos<double, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

BOOST_AUTO_TEST_CASE( gamma_ele_pos_6 )
{
    const double px = -2534.83*me_c;
    const double py = 0*me_c;
    const double pz = 0*me_c;

    const double gamma_exp = 2534.830197251879;

    const auto p = vec3<double>{px,py,pz};

    test_gamma_ele_pos<double, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<double, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<double, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<double, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::SI>(p, gamma_exp, 1.0);
    test_gamma_ele_pos<float, unit_system::norm_lambda>(p, gamma_exp, lambda);
    test_gamma_ele_pos<float, unit_system::norm_omega>(p, gamma_exp, omega);
    test_gamma_ele_pos<float, unit_system::heaviside_lorentz>(p, gamma_exp, 1.0);
}

//*******************************
