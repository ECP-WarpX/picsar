//####### Test module for Quantum Synchrotron emission engine (SI units) #######

//Define Module name
 #define BOOST_TEST_MODULE "Quantum Synchrotron emission engine (SI units)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include <utility>

 //Include Boost unit tests library & library for floating point comparison
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>

//SI units are used for this test
#define PXRMP_USE_SI_UNITS
#include "quantum_sync_engine.hpp"

using namespace picsar::multi_physics;

//________________________________

//Helper function
template<typename REAL>
quantum_synchrotron_engine<REAL, stl_rng_wrapper> get_qs_stl_set_lambda(uint64_t seed, REAL lambda,
quantum_synchrotron_engine_ctrl<REAL> qs_ctrl = quantum_synchrotron_engine_ctrl<REAL>())
{
    stl_rng_wrapper wrap{seed};
    auto qs_engine =  quantum_synchrotron_engine<REAL, stl_rng_wrapper>{std::move(wrap), 1.0, qs_ctrl};
    qs_engine.set_lambda(static_cast<REAL>(lambda));
    return qs_engine;
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


//Test get/set lambda for quantum_synchrotron_engine generic
template <typename T>
void quantum_sync_engine_gs()
{
    auto qs_engine = get_qs_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));
    BOOST_CHECK_EQUAL( static_cast<T>(1.0), qs_engine.get_lambda());
}

//Test get/set lambda for quantum_synchrotron_engine (double precision)
BOOST_AUTO_TEST_CASE( quantum_sync_engine_gs_double_1 )
{
    quantum_sync_engine_gs<double>();
}

//Test get/set lambda for quantum_synchrotron_engine (single precision)
BOOST_AUTO_TEST_CASE( quantum_sync_engine_gs_single_1 )
{
    quantum_sync_engine_gs<float>();
}

// ------------- optical depth --------------
//Test get/set lambda for quantum_synchrotron_engine generic
template <typename T>
void quantum_sync_engine_opt()
{
    auto qs_engine = get_qs_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));
    BOOST_TEST ( qs_engine.get_optical_depth() >= static_cast<T>(0.0) );
}

//Test get new optical depth (double precision)
BOOST_AUTO_TEST_CASE( quantum_sync_engine_opt_double_1 )
{
    quantum_sync_engine_opt<double>();
}

//Test get new optical depth (single precision)
BOOST_AUTO_TEST_CASE( quantum_sync_engine_opt_single_1 )
{
    quantum_sync_engine_opt<float>();
}
