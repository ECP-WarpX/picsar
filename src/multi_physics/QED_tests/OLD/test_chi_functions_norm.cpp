//####### Test module for chi functions (normalized units) #####################

//Define Module name
 #define BOOST_TEST_MODULE "chi functions (normalized units)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//For this test we will use normalized units
#define PXRMP_USE_NORMALIZED_UNITS
#include "chi_functions.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-10;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//#################### Photons

//Test chi function for photons with double precision (case 1)
BOOST_AUTO_TEST_CASE( chi_photons_double_1 )
{
    double px = 83.759;
    double py = 139.311;
    double pz = -230.553;
    double ex = -166.145;
    double ey = -78.231;
    double ez = -278.856;
    double bx = -279.174;
    double by = -158.849;
    double bz = -93.826;
    double lambda = 800. * si_nanometer;

    double chi_exp = 0.347111844317;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for photons with double precision (case 2)
BOOST_AUTO_TEST_CASE( chi_photons_double_2 )
{
    double px = 9.2627;
    double py = -25.4575;
    double pz = -10.2246;
    double ex = 2.9271;
    double ey = 10.4293;
    double ez = 3.6103;
    double bx = 1.7439;
    double by = 1.9778;
    double bz = 17.8799;
    double lambda = 800. * si_nanometer;

    double chi_exp = 0.00090414740533;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}


//Test chi function for photons with double precision (case 3)
BOOST_AUTO_TEST_CASE( chi_photons_double_3 )
{
    double px = -2314.45;
    double py = -2356.30;
    double pz = 546.28;
    double ex = 1230.11;
    double ey =  1638.02;
    double ez = -2911.04;
    double bx = -2203.66;
    double by = 1243.79;
    double bz = -2830.99;
    double lambda = 800. * si_nanometer;

    double chi_exp = 57.2204397969;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for photons with double precision (case 4)
BOOST_AUTO_TEST_CASE( chi_photons_double_4 )
{
    double px = 0;
    double py = 0;
    double pz = 0;
    double ex = 1230.11;
    double ey =  1638.02;
    double ez = -2911.04;
    double bx = -2203.66;
    double by = 1243.79;
    double bz = -2830.99;
    double lambda = 800. * si_nanometer;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for photons with double precision (case 5)
BOOST_AUTO_TEST_CASE( chi_photons_double_5 )
{
    double px = -2314.45;
    double py = -2356.30;
    double pz = 546.28;
    double ex = 0;
    double ey =  0;
    double ez = 0;
    double bx = 0;
    double by = 0;
    double bz = 0;
    double lambda = 800. * si_nanometer;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for photons with double precision (case 6)
BOOST_AUTO_TEST_CASE( chi_photons_double_6 )
{
    double px = -2314.45;
    double py = 0;
    double pz = 0;
    double ex = 1230.11;
    double ey = 0;
    double ez = 0;
    double bx = 0;
    double by = 0;
    double bz = 0;
    double lambda = 800. * si_nanometer;

    double chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL(chi_res, double_tolerance);
}

//Test chi function for photons with single precision (case 1)
BOOST_AUTO_TEST_CASE( chi_photons_single_1 )
{
    float px = 83.759f;
    float py = 139.311f;
    float pz = -230.553f;
    float ex = -166.145f;
    float ey = -78.231f;
    float ez = -278.856f;
    float bx = -279.174f;
    float by = -158.849f;
    float bz = -93.826f;
    float lambda = 800.0f * si_nanometer;

    float chi_exp = static_cast<float>(0.347111844317);

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for photons with single precision (case 2)
BOOST_AUTO_TEST_CASE( chi_photons_single_2 )
{
    float px = 9.2627f;
    float py = -25.4575f;
    float pz = -10.2246f;
    float ex = 2.9271f;
    float ey = 10.4293f;
    float ez = 3.6103f;
    float bx = 1.7439f;
    float by = 1.9778f;
    float bz = 17.8799f;
    float lambda = 800.0f * si_nanometer;

    float chi_exp = static_cast<float>(0.00090414740533);

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}


//Test chi function for photons with single precision (case 3)
BOOST_AUTO_TEST_CASE( chi_photons_single_3 )
{
    float px = -2314.45f;
    float py = -2356.30f;
    float pz = 546.28f;
    float ex = 1230.11f;
    float ey =  1638.02f;
    float ez = -2911.04f;
    float bx = -2203.66f;
    float by = 1243.79f;
    float bz = -2830.99f;
    float lambda = 800.0f * si_nanometer;

    float chi_exp = static_cast<float>(57.2204397969);

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for photons with single precision (case 4)
BOOST_AUTO_TEST_CASE( chi_photons_float_4 )
{
    float px = 0.0f;
    float py = 0.0f;
    float pz = 0.0f;
    float ex = 1230.11f;
    float ey =  1638.02f;
    float ez = -2911.04f;
    float bx = -2203.66f;
    float by = 1243.79f;
    float bz = -2830.99f;
    float lambda = 800.0f * si_nanometer;

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for photons with single precision (case 5)
BOOST_AUTO_TEST_CASE( chi_photons_float_5 )
{
    float px = -2314.45f;
    float py = -2356.30f;
    float pz = 546.28f;
    float ex = 0.0f;
    float ey =  0.0f;
    float ez = 0.0f;
    float bx = 0.0f;
    float by = 0.0f;
    float bz = 0.0f;
    float lambda = 800.0f * si_nanometer;

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for photons with single precision (case 6)
BOOST_AUTO_TEST_CASE( chi_photons_float_6 )
{
    float px = -2314.45f;
    float py = 0.0f;
    float pz = 0.0f;
    float ex = 1230.11f;
    float ey = 0.0f;
    float ez = 0.0f;
    float bx = 0.0f;
    float by = 0.0f;
    float bz = 0.0f;
    float lambda = 800.0f * si_nanometer;

    float chi_res = chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL(chi_res, float_tolerance);
}

//#################### Leptons

//Test chi function for leptons with double precision (case 1)
BOOST_AUTO_TEST_CASE( chi_leptons_double_1 )
{
    double px = 24.3752;
    double py = -11.5710;
    double pz = -10.0841;
    double ex = 57.185;
    double ey = -16.6555;
    double ez = 22.4340;
    double bx = 6.6911;
    double by = -23.8724;
    double bz = 13.9934;
    double lambda = 800. * si_nanometer;

    double chi_exp = 0.00216716627219670;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for leptons with double precision (case 2)
BOOST_AUTO_TEST_CASE( chi_leptons_double_2 )
{
    double px = 4.015;
    double py = 197.287;
    double pz = 141.705;
    double ex = 30.287;
    double ey = 115.740;
    double ez = 120.891;
    double bx = -190.161;
    double by = -129.115;
    double bz = -57.002;
    double lambda = 800. * si_nanometer;

    double chi_exp = 0.166318112874468;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for leptons with double precision (case 3)
BOOST_AUTO_TEST_CASE( chi_leptons_double_3 )
{
    double px = -2534.83;
    double py = 1011.54;
    double pz = -793.04;
    double ex = 741.67;
    double ey = -2359.97;
    double ez = 1463.50;
    double bx = 1477.19;
    double by = -1448.33;
    double bz = 1953.68;
    double lambda = 800. * si_nanometer;

    double chi_exp = 16.0114572646993;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, double_tolerance);
}

//Test chi function for leptons with double precision (case 4)
BOOST_AUTO_TEST_CASE( chi_leptons_double_4 )
{
    double px = 0;
    double py = 0;
    double pz = 0;
    double ex = 741.67;
    double ey = -2359.97;
    double ez = 1463.50;
    double bx = 1477.19;
    double by = -1448.33;
    double bz = 1953.68;
    double lambda = 800. * si_nanometer;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for leptons with double precision (case 5)
BOOST_AUTO_TEST_CASE( chi_leptons_double_5 )
{
    double px = -2534.83;
    double py = 1011.54;
    double pz = -793.04;
    double ex = 0;
    double ey =  0;
    double ez = 0;
    double bx = 0;
    double by = 0;
    double bz = 0;
    double lambda = 800. * si_nanometer;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0);
}

//Test chi function for leptons with double precision (case 6)
BOOST_AUTO_TEST_CASE( chi_leptons_double_6 )
{
    double px = -2534.83;
    double py = 0;
    double pz = 0;
    double ex = 0;
    double ey = 0;
    double ez = 0;
    double bx = 1477.19;
    double by = 0;
    double bz = 0;
    double lambda = 800. * si_nanometer;

    double chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL(chi_res, double_tolerance);
}

//Test chi function for leptons with single precision (case 1)
BOOST_AUTO_TEST_CASE( chi_leptons_single_1 )
{
    float px = 24.3752f;
    float py = -11.5710f;
    float pz = -10.0841f;
    float ex = 57.185f;
    float ey = -16.6555f;
    float ez = 22.4340f;
    float bx = 6.6911f;
    float by = -23.8724f;
    float bz = 13.9934f;
    float lambda = 800.0f * si_nanometer;

    float chi_exp = static_cast<float>(0.00216716627219670);

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for leptons with single precision (case 2)
BOOST_AUTO_TEST_CASE( chi_leptons_single_2 )
{
    float px = 4.015f;
    float py = 197.287f;
    float pz = 141.705f;
    float ex = 30.287f;
    float ey = 115.740f;
    float ez = 120.891f;
    float bx = -190.161f;
    float by = -129.115f;
    float bz = -57.002f;
    float lambda = 800.0f * si_nanometer;

    float chi_exp = static_cast<float>(0.166318112874468);

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}


//Test chi function for leptons with single precision (case 3)
BOOST_AUTO_TEST_CASE( chi_leptons_single_3 )
{
    float px = -2534.83f;
    float py = 1011.54f;
    float pz = -793.04f;
    float ex = 741.67f;
    float ey = -2359.97f;
    float ez = 1463.50f;
    float bx = 1477.19f;
    float by = -1448.33f;
    float bz = 1953.68f;
    float lambda = 800.0f * si_nanometer;

    float chi_exp = static_cast<float>(16.0114572646993);

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL((chi_res-chi_exp)/chi_exp, float_tolerance);
}

//Test chi function for leptons with single precision (case 4)
BOOST_AUTO_TEST_CASE( chi_leptons_single_4 )
{
    float px = 0.0f;
    float py = 0.0f;
    float pz = 0.0f;
    float ex = 741.67f;
    float ey = -2359.97f;
    float ez = 1463.50f;
    float bx = 1477.19f;
    float by = -1448.33f;
    float bz = 1953.68f;
    float lambda = 800.0f * si_nanometer;

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for leptons with single precision (case 5)
BOOST_AUTO_TEST_CASE( chi_leptons_single_5 )
{
    float px = -2534.83f;
    float py = 1011.54f;
    float pz = -793.04f;
    float ex = 0.0f;
    float ey = 0.0f;
    float ez = 0.0f;
    float bx = 0.0f;
    float by = 0.0f;
    float bz = 0.0f;
    float lambda = 800.0f * si_nanometer;

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_EQUAL(chi_res, 0.0f);
}

//Test chi function for leptons with single precision (case 6)
BOOST_AUTO_TEST_CASE( chi_leptons_single_6 )
{
    float px = -2534.83f;
    float py = 0.0f;
    float pz = 0.0f;
    float ex = 0.0f;
    float ey = 0.0f;
    float ez = 0.0f;
    float bx = 1477.19f;
    float by = 0.0f;
    float bz = 0.0f;
    float lambda = 800.0f * si_nanometer;

    float chi_res = chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda);

    BOOST_CHECK_SMALL(chi_res, float_tolerance);
}
