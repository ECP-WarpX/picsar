//####### Test module for chi functions (SI units) #############################

//Define Module name
 #define BOOST_TEST_MODULE "chi functions (SI units)"

//Will automatically define a main for this test
#define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//For this test we will use SI units
#define PXRMP_USE_SI_UNITS
#include "chi_functions.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-10;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//SI units for momenta
const double me_c = electron_mass * light_speed;
const float flt_me_c = flt_electron_mass * flt_light_speed;

//SI units for fields
double lambda = 800.0 * si_nanometer;
double eref = 2.0*pi*electron_mass*light_speed*light_speed/
            (lambda*elementary_charge);
double bref = eref/light_speed;

float flt_lambda = 800.0f * flt_si_nanometer;
float flt_eref = 2.0f*flt_pi*flt_electron_mass*flt_light_speed*flt_light_speed/
            (flt_lambda*flt_elementary_charge);
float flt_bref = flt_eref/flt_light_speed;

//#################### Photons

//Test chi function for photons with double precision (case 1)
BOOST_AUTO_TEST_CASE( chi_photons_double_1 )
{
    double px = 83.759*me_c;
    double py = 139.311*me_c;
    double pz = -230.553*me_c;
    double ex = -166.145*eref;
    double ey = -78.231*eref;
    double ez = -278.856*eref;
    double bx = -279.174*bref;
    double by = -158.849*bref;
    double bz = -93.826*bref;

    double chi_exp = 0.347111844317;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for photons with double precision (case 2)
BOOST_AUTO_TEST_CASE( chi_photons_double_2 )
{
    double px = 9.2627*me_c;
    double py = -25.4575*me_c;
    double pz = -10.2246*me_c;
    double ex = 2.9271*eref;
    double ey = 10.4293*eref;
    double ez = 3.6103*eref;
    double bx = 1.7439*bref;
    double by = 1.9778*bref;
    double bz = 17.8799*bref;

    double chi_exp = 0.00090414740533;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}


//Test chi function for photons with double precision (case 3)
BOOST_AUTO_TEST_CASE( chi_photons_double_3 )
{
    double px = -2314.45*me_c;
    double py = -2356.30*me_c;
    double pz = 546.28*me_c;
    double ex = 1230.11*eref;
    double ey =  1638.02*eref;
    double ez = -2911.04*eref;
    double bx = -2203.66*bref;
    double by = 1243.79*bref;
    double bz = -2830.99*bref;

    double chi_exp = 57.2204397969;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for photons with double precision (case 4)
BOOST_AUTO_TEST_CASE( chi_photons_double_4 )
{
    double px = 0*me_c;
    double py = 0*me_c;
    double pz = 0*me_c;
    double ex = 1230.11*eref;
    double ey =  1638.02*eref;
    double ez = -2911.04*eref;
    double bx = -2203.66*bref;
    double by = 1243.79*bref;
    double bz = -2830.99*bref;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for photons with double precision (case 5)
BOOST_AUTO_TEST_CASE( chi_photons_double_5 )
{
    double px = -2314.45*me_c;
    double py = -2356.30*me_c;
    double pz = 546.28*me_c;
    double ex = 0*eref;
    double ey =  0*eref;
    double ez = 0*eref;
    double bx = 0*bref;
    double by = 0*bref;
    double bz = 0*bref;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for photons with double precision (case 6)
BOOST_AUTO_TEST_CASE( chi_photons_double_6 )
{
    double px = -2314.45*me_c;
    double py = 0*me_c;
    double pz = 0*me_c;
    double ex = 1230.11*eref;
    double ey = 0*eref;
    double ez = 0*eref;
    double bx = 0*bref;
    double by = 0*bref;
    double bz = 0*bref;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL(chi_res, double_tolerance);
}

//Test chi function for photons with single precision (case 1)
BOOST_AUTO_TEST_CASE( chi_photons_single_1 )
{
    float px = 83.759*flt_me_c;
    float py = 139.311f*flt_me_c;
    float pz = -230.553f*flt_me_c;
    float ex = -166.145f*flt_eref;
    float ey = -78.231f*flt_eref;
    float ez = -278.856f*flt_eref;
    float bx = -279.174f*flt_bref;
    float by = -158.849f*flt_bref;
    float bz = -93.826f*flt_bref;

    float chi_exp = static_cast<float>(0.347111844317);

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for photons with single precision (case 2)
BOOST_AUTO_TEST_CASE( chi_photons_single_2 )
{
    float px = 9.2627f*flt_me_c;
    float py = -25.4575f*flt_me_c;
    float pz = -10.2246f*flt_me_c;
    float ex = 2.9271f*flt_eref;
    float ey = 10.4293f*flt_eref;
    float ez = 3.6103f*flt_eref;
    float bx = 1.7439f*flt_bref;
    float by = 1.9778f*flt_bref;
    float bz = 17.8799f*flt_bref;

    float chi_exp = static_cast<float>(0.00090414740533);

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}


//Test chi function for photons with single precision (case 3)
BOOST_AUTO_TEST_CASE( chi_photons_single_3 )
{
    float px = -2314.45f*flt_me_c;
    float py = -2356.30f*flt_me_c;
    float pz = 546.28f*flt_me_c;
    float ex = 1230.11f*flt_eref;
    float ey =  1638.02f*flt_eref;
    float ez = -2911.04f*flt_eref;
    float bx = -2203.66f*flt_bref;
    float by = 1243.79f*flt_bref;
    float bz = -2830.99f*flt_bref;

    float chi_exp = static_cast<float>(57.2204397969);

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for photons with single precision (case 4)
BOOST_AUTO_TEST_CASE( chi_photons_float_4 )
{
    float px = 0.0f*flt_me_c;
    float py = 0.0f*flt_me_c;
    float pz = 0.0f*flt_me_c;
    float ex = 1230.11f*flt_eref;
    float ey =  1638.02f*flt_eref;
    float ez = -2911.04f*flt_eref;
    float bx = -2203.66f*flt_bref;
    float by = 1243.79f*flt_bref;
    float bz = -2830.99f*flt_bref;

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for photons with single precision (case 5)
BOOST_AUTO_TEST_CASE( chi_photons_float_5 )
{
    float px = -2314.45f*flt_me_c;
    float py = -2356.30f*flt_me_c;
    float pz = 546.28f*flt_me_c;
    float ex = 0.0f*flt_eref;
    float ey =  0.0f*flt_eref;
    float ez = 0.0f*flt_eref;
    float bx = 0.0f*flt_bref;
    float by = 0.0f*flt_bref;
    float bz = 0.0f*flt_bref;

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for photons with single precision (case 6)
BOOST_AUTO_TEST_CASE( chi_photons_float_6 )
{
    float px = -2314.45f*flt_me_c;
    float py = 0.0f*flt_me_c;
    float pz = 0.0f*flt_me_c;
    float ex = 1230.11f*flt_eref;
    float ey = 0.0f*flt_eref;
    float ez = 0.0f*flt_eref;
    float bx = 0.0f*flt_bref;
    float by = 0.0f*flt_bref;
    float bz = 0.0f*flt_bref;

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL(chi_res, float_tolerance);
}

//#################### Leptons

//Test chi function for leptons with double precision (case 1)
BOOST_AUTO_TEST_CASE( chi_leptons_double_1 )
{
    double px = 24.3752*me_c;
    double py = -11.5710*me_c;
    double pz = -10.0841*me_c;
    double ex = 57.185*eref;
    double ey = -16.6555*eref;
    double ez = 22.4340*eref;
    double bx = 6.6911*bref;
    double by = -23.8724*bref;
    double bz = 13.9934*bref;

    double chi_exp = 0.00216716627219670;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for leptons with double precision (case 2)
BOOST_AUTO_TEST_CASE( chi_leptons_double_2 )
{
    double px = 4.015*me_c;
    double py = 197.287*me_c;
    double pz = 141.705*me_c;
    double ex = 30.287*eref;
    double ey = 115.740*eref;
    double ez = 120.891*eref;
    double bx = -190.161*bref;
    double by = -129.115*bref;
    double bz = -57.002*bref;

    double chi_exp = 0.166318112874468;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for leptons with double precision (case 3)
BOOST_AUTO_TEST_CASE( chi_leptons_double_3 )
{
    double px = -2534.83*me_c;
    double py = 1011.54*me_c;
    double pz = -793.04*me_c;
    double ex = 741.67*eref;
    double ey = -2359.97*eref;
    double ez = 1463.50*eref;
    double bx = 1477.19*bref;
    double by = -1448.33*bref;
    double bz = 1953.68*bref;

    double chi_exp = 16.0114572646993;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for leptons with double precision (case 4)
BOOST_AUTO_TEST_CASE( chi_leptons_double_4 )
{
    double px = 0*me_c;
    double py = 0*me_c;
    double pz = 0*me_c;
    double ex = 741.67*eref;
    double ey = -2359.97*eref;
    double ez = 1463.50*eref;
    double bx = 1477.19*bref;
    double by = -1448.33*bref;
    double bz = 1953.68*bref;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for leptons with double precision (case 5)
BOOST_AUTO_TEST_CASE( chi_leptons_double_5 )
{
    double px = -2534.83*me_c;
    double py = 1011.54*me_c;
    double pz = -793.04*me_c;
    double ex = 0*eref;
    double ey = 0*eref;
    double ez = 0*eref;
    double bx = 0*bref;
    double by = 0*bref;
    double bz = 0*bref;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for leptons with double precision (case 6)
BOOST_AUTO_TEST_CASE( chi_leptons_double_6 )
{
    double px = -2534.83*me_c;
    double py = 0*me_c;
    double pz = 0*me_c;
    double ex = 0*eref;
    double ey = 0*eref;
    double ez = 0*eref;
    double bx = 1477.19*bref;
    double by = 0*bref;
    double bz = 0*bref;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL(chi_res, double_tolerance);
}

//Test chi function for leptons with single precision (case 1)
BOOST_AUTO_TEST_CASE( chi_leptons_single_1 )
{
    float px = 24.3752f*flt_me_c;
    float py = -11.5710f*flt_me_c;
    float pz = -10.0841f*flt_me_c;
    float ex = 57.185f*flt_eref;
    float ey = -16.6555f*flt_eref;
    float ez = 22.4340f*flt_eref;
    float bx = 6.6911f*flt_bref;
    float by = -23.8724f*flt_bref;
    float bz = 13.9934f*flt_bref;

    float chi_exp = static_cast<float>(0.00216716627219670);

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for leptons with single precision (case 2)
BOOST_AUTO_TEST_CASE( chi_leptons_single_2 )
{
    float px = 4.015f*flt_me_c;
    float py = 197.287f*flt_me_c;
    float pz = 141.705f*flt_me_c;
    float ex = 30.287f*flt_eref;
    float ey = 115.740f*flt_eref;
    float ez = 120.891f*flt_eref;
    float bx = -190.161f*flt_bref;
    float by = -129.115f*flt_bref;
    float bz = -57.002f*flt_bref;

    float chi_exp = static_cast<float>(0.166318112874468);

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}


//Test chi function for leptons with single precision (case 3)
BOOST_AUTO_TEST_CASE( chi_leptons_single_3 )
{
    float px = -2534.83f*flt_me_c;
    float py = 1011.54f*flt_me_c;
    float pz = -793.04f*flt_me_c;
    float ex = 741.67f*flt_eref;
    float ey = -2359.97f*flt_eref;
    float ez = 1463.50f*flt_eref;
    float bx = 1477.19f*flt_bref;
    float by = -1448.33f*flt_bref;
    float bz = 1953.68f*flt_bref;

    float chi_exp = static_cast<float>(16.0114572646993);

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for leptons with single precision (case 4)
BOOST_AUTO_TEST_CASE( chi_leptons_single_4 )
{
    float px = 0.0f*flt_me_c;
    float py = 0.0f*flt_me_c;
    float pz = 0.0f*flt_me_c;
    float ex = 741.67f*flt_eref;
    float ey = -2359.97f*flt_eref;
    float ez = 1463.50f*flt_eref;
    float bx = 1477.19f*flt_bref;
    float by = -1448.33f*flt_bref;
    float bz = 1953.68f*flt_bref;

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for leptons with single precision (case 5)
BOOST_AUTO_TEST_CASE( chi_leptons_single_5 )
{
    float px = -2534.83f*flt_me_c;
    float py = 1011.54f*flt_me_c;
    float pz = -793.04f*flt_me_c;
    float ex = 0.0f*flt_eref;
    float ey = 0.0f*flt_eref;
    float ez = 0.0f*flt_eref;
    float bx = 0.0f*flt_bref;
    float by = 0.0f*flt_bref;
    float bz = 0.0f*flt_bref;

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for leptons with single precision (case 6)
BOOST_AUTO_TEST_CASE( chi_leptons_single_6 )
{
    float px = -2534.83f*flt_me_c;
    float py = 0.0f*flt_me_c;
    float pz = 0.0f*flt_me_c;
    float ex = 0.0f*flt_eref;
    float ey = 0.0f*flt_eref;
    float ez = 0.0f*flt_eref;
    float bx = 1477.19f*flt_bref;
    float by = 0.0f*flt_bref;
    float bz = 0.0f*flt_bref;

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz);

    BOOST_CHECK_SMALL(chi_res, float_tolerance);
}
