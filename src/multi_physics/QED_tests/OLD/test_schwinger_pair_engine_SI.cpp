//####### Test module for Schwinger Pair engine (SI units) ######################

//Define Module name
 #define BOOST_TEST_MODULE "Schwinger Pair engine engine (SI)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include <memory>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//SI units are used for this test
#define PXRMP_USE_SI_UNITS
#include "schwinger_pair_engine.hpp"

#include "rng_wrapper.hpp"

using namespace picsar::multi_physics;

//Helper function
template<typename REAL>
schwinger_pair_engine<REAL, stl_rng_wrapper<REAL>>
get_schwinger_stl_set_lambda(uint64_t seed)
{
    stl_rng_wrapper<REAL> wrap{seed};
    auto schw_engine =  schwinger_pair_engine<REAL, stl_rng_wrapper<REAL>>
        {std::move(wrap), 1.0};
    return schw_engine;
}

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-3;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-2;

//Templated tolerance
template <typename T>
T tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//***SI UNITS***
//SI units for momenta
const double me_c = electron_mass * light_speed;
//SI units for fields
double lambda = 800.0 * si_nanometer;
double eref = 2.0*pi*electron_mass*light_speed*light_speed/
            (lambda*elementary_charge);
double bref = eref/light_speed;
//SI units for dt and rate
double dtref = lambda/(2.0*pi*light_speed);
double rateref = 1.0/dtref;

//Test get/set lambda for breit_wheeler_engine generic
template <typename T>
void schwinger_pair_engine_gs()
{
    auto schw_engine = get_schwinger_stl_set_lambda<T>(390109317);
    BOOST_CHECK_EQUAL( static_cast<T>(1.0), schw_engine.get_lambda());
}

//Test get/set lambda for breit_wheeler_engine (double precision)
BOOST_AUTO_TEST_CASE( schwinger_pair_engine_gs_double_1 )
{
    schwinger_pair_engine_gs<double>();
}

//Test get/set lambda for breit_wheeler_engine (single precision)
BOOST_AUTO_TEST_CASE( schwinger_pair_engine_gs_single_1 )
{
    schwinger_pair_engine_gs<float>();
}

// ------------- Pair positions --------------
//Test relative position of the new pair for the schwinger_pair_engine
template <typename T>
void schwinger_pair_engine_pos()
{
    auto schw_engine = get_schwinger_stl_set_lambda<T>(390109317);
    auto res = schw_engine.get_new_pair_position();
    for (auto val: res)
        BOOST_TEST ( (val >= static_cast<T>(0.0) && val <= static_cast<T>(1.0)) );
}

//Test get new optical depth (double precision)
BOOST_AUTO_TEST_CASE( schwinger_pair_engine_pos_double_1 )
{
    schwinger_pair_engine_pos<double>();
}

//Test get new optical depth (single precision)
BOOST_AUTO_TEST_CASE( schwinger_pair_engine_pos_single_1 )
{
    schwinger_pair_engine_pos<float>();
}

// ------------- Pair production rates --------------

/*
//Test pair production rates (generic)
template <typename T>
void breit_wheeler_engine_prod_6()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(-965.61*me_c);
    T py =  static_cast<T>(-3975.11*me_c);
    T pz =  static_cast<T>(6917.22*me_c);
    T ex =  static_cast<T>(11.17*eref);
    T ey =  static_cast<T>(-2117.72*eref);
    T ez =  static_cast<T>(-1407.19*eref);
    T bx =  static_cast<T>( 6259.79*bref);
    T by =  static_cast<T>(7557.54*bref);
    T bz =  static_cast<T>(773.11*bref);

    T exp = static_cast<T>(3.51855878777*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL((exp-res)/exp, tolerance<T>());
}

//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_6 )
{
    breit_wheeler_engine_prod_6<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_6 )
{
    breit_wheeler_engine_prod_6<float>();
}
*/
