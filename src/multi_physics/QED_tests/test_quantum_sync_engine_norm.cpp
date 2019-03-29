//####### Test module for Quantum Synchrotron emission engine (normalized units) #######

//Define Module name
 #define BOOST_TEST_MODULE "Quantum Synchrotron emission engine (normalized units)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include <utility>

 //Include Boost unit tests library & library for floating point comparison
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>

//Normalized units are used for this test
#define PXRMP_USE_NORMALIZED_UNITS
#include "quantum_sync_engine.hpp"


//Only include Kokkos if it is supported
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT
    #include <Kokkos_Core.hpp>
    #include <Kokkos_Random.hpp>
    #include <Kokkos_DualView.hpp>
#endif

using namespace picsar::multi_physics;

//________________________________

//Helper function
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

//Kokkos version
template<typename REAL>
quantum_synchrotron_engine
<REAL, kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>>>
get_qs_kokkos_set_lambda(uint64_t seed, REAL lambda,
quantum_synchrotron_engine_ctrl<REAL> qs_ctrl = quantum_synchrotron_engine_ctrl<REAL>())
{
    auto pool = Kokkos::Random_XorShift1024_Pool<>{seed};
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrap{pool};
    auto qs_engine =  quantum_synchrotron_engine
        <REAL, kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>>>{std::move(wrap), 1.0, qs_ctrl};
    qs_engine.set_lambda(static_cast<REAL>(lambda));
    return qs_engine;
}

#endif

//No Kokkos version
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

//Test get/set lambda for quantum_synchrotron_engine generic
template <typename T>
void quantum_sync_engine_gs()
{
    auto qs_engine = get_qs_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));
    BOOST_CHECK_EQUAL( static_cast<T>(800.0*si_nanometer), qs_engine.get_lambda());
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
//Test get/set lambda forquantum_synchrotron_engine generic
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

// ------------- photon emission rate --------------

//Test photon emission (generic)
template <typename T>
void quantum_sync_engine_phtoem()
{
    T lambda = static_cast<T>(800. * si_nanometer);

    auto qs_engine = get_qs_stl_set_lambda<T>(390109317,lambda);

    T px =  static_cast<T>(149.825);
    T py =  static_cast<T>(933.115);
    T pz =  static_cast<T>(-538.195);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(-861.074);
    T ez =  static_cast<T>(944.652);
    T bx =  static_cast<T>(531.406);
    T by =  static_cast<T>(670.933);
    T bz =  static_cast<T>(660.057);

    //T exp = static_cast<T>(xxx);
    T res = qs_engine.compute_dN_dt
        (chi_lepton(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );


    BOOST_CHECK_EQUAL(1,1);
    //BOOST_CHECK_SMALL(res, tolerance<T>());
}

//Test photon emission (double precision)
BOOST_AUTO_TEST_CASE( quantum_sync_engine_phtoem_double_1 )
{
    quantum_sync_engine_phtoem<double>();
}

//Test photon emission (single precision)
BOOST_AUTO_TEST_CASE( quantum_sync_engine_phtoem_single_1 )
{
    quantum_sync_engine_phtoem<float>();
}
