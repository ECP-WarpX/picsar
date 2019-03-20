//####### Test module for Breit Wheeler engine (SI units) ######################

//Define Module name
 #define BOOST_TEST_MODULE "Nonliner Breit Wheeler engine (SI)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include <memory>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//SI units are used for this test
#define PXRMP_USE_SI_UNITS
#include "breit_wheeler_engine.hpp"

#include "chi_functions.hpp"
#include "vec_functions.hpp"

#include "rng_wrapper.hpp"

using namespace picsar::multi_physics;

//Helper function
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

//Kokkos version
template<typename REAL>
breit_wheeler_engine
<REAL, kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>>>
get_bw_kokkos_set_lambda(uint64_t seed,
breit_wheeler_engine_ctrl<REAL> bw_ctrl = breit_wheeler_engine_ctrl<REAL>())
{
    auto pool = Kokkos::Random_XorShift1024_Pool<>{seed};
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrap{pool};
    auto bw_engine =  breit_wheeler_engine
        <REAL, kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>>>{std::move(wrap), 1.0, bw_ctrl};
    return bw_engine;
}

#endif

//No Kokkos version
template<typename REAL>
breit_wheeler_engine<REAL, stl_rng_wrapper> get_bw_stl_set_lambda(uint64_t seed,
breit_wheeler_engine_ctrl<REAL> bw_ctrl = breit_wheeler_engine_ctrl<REAL>())
{
    stl_rng_wrapper wrap{seed};
    auto bw_engine =  breit_wheeler_engine<REAL, stl_rng_wrapper>{std::move(wrap), 1.0, bw_ctrl};
    return bw_engine;
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
void breit_wheeler_engine_gs()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);
    BOOST_CHECK_EQUAL( static_cast<T>(1.0), bw_engine.get_lambda());
}

//Test get/set lambda for breit_wheeler_engine (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_gs_double_1 )
{
    breit_wheeler_engine_gs<double>();
}

//Test get/set lambda for breit_wheeler_engine (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_gs_single_1 )
{
    breit_wheeler_engine_gs<float>();
}

// ------------- optical depth --------------
//Test get/set lambda for breit_wheeler_engine generic
template <typename T>
void breit_wheeler_engine_opt()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);
    BOOST_TEST ( bw_engine.get_optical_depth() >= static_cast<T>(0.0) );
}

//Test get new optical depth (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_opt_double_1 )
{
    breit_wheeler_engine_opt<double>();
}

//Test get new optical depth (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_opt_single_1 )
{
    breit_wheeler_engine_opt<float>();
}

// ------------- Pair production rates --------------

//Test pair production rates (generic)
template <typename T>
void breit_wheeler_engine_prod_1()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(0.0*me_c);
    T py =  static_cast<T>(0.0*me_c);
    T pz =  static_cast<T>(0.0*me_c);
    T ex =  static_cast<T>(931.686*eref);
    T ey =  static_cast<T>(-861.074*eref);
    T ez =  static_cast<T>(944.652*eref);
    T bx =  static_cast<T>(531.406*bref);
    T by =  static_cast<T>(670.933*bref);
    T bz =  static_cast<T>(660.057*bref);

    //T exp = static_cast<T>(0.0*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL(res, tolerance<T>());
}


//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_1 )
{
    breit_wheeler_engine_prod_1<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_1 )
{
    breit_wheeler_engine_prod_1<float>();
}

//Test pair production rates (generic)
template <typename T>
void breit_wheeler_engine_prod_2()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(149.825*me_c);
    T py =  static_cast<T>(0.0*me_c);
    T pz =  static_cast<T>(0.0*me_c);
    T ex =  static_cast<T>(931.686*eref);
    T ey =  static_cast<T>(0.0*eref);
    T ez =  static_cast<T>(0.0*eref);
    T bx =  static_cast<T>(0.0*bref);
    T by =  static_cast<T>(0.0*bref);
    T bz =  static_cast<T>(0.0*bref);

    //T exp = static_cast<T>(0.0*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL(res, tolerance<T>());
}

//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_2 )
{
    breit_wheeler_engine_prod_2<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_2 )
{
    breit_wheeler_engine_prod_2<float>();
}

template <typename T>
void breit_wheeler_engine_prod_3()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(6.81696*me_c);
    T py =  static_cast<T>(9.68933*me_c);
    T pz =  static_cast<T>(2.81229*me_c);
    T ex =  static_cast<T>(-4.89986*eref);
    T ey =  static_cast<T>(-9.65535*eref);
    T ez =  static_cast<T>(3.69471*eref);
    T bx =  static_cast<T>(8.89549*bref);
    T by =  static_cast<T>(-5.46574*bref);
    T bz =  static_cast<T>(-6.75393*bref);

    //T exp = static_cast<T>(0.0*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL(res, tolerance<T>());
}


//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_3 )
{
    breit_wheeler_engine_prod_3<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_3 )
{
    breit_wheeler_engine_prod_3<float>();
}

//Test pair production rates (generic)
template <typename T>
void breit_wheeler_engine_prod_4()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(149.825*me_c);
    T py =  static_cast<T>(933.115*me_c);
    T pz =  static_cast<T>(-538.195*me_c);
    T ex =  static_cast<T>(931.686*eref);
    T ey =  static_cast<T>(-861.074*eref);
    T ez =  static_cast<T>(944.652*eref);
    T bx =  static_cast<T>(531.406*bref);
    T by =  static_cast<T>(670.933*bref);
    T bz =  static_cast<T>(660.057*bref);

    T exp = static_cast<T>(1.50648551484*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL((exp-res)/exp, tolerance<T>());
}


//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_4 )
{
    breit_wheeler_engine_prod_4<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_4 )
{
    breit_wheeler_engine_prod_4<float>();
}



//Test pair production rates (generic)
template <typename T>
void breit_wheeler_engine_prod_5()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(-44.4546*me_c);
    T py =  static_cast<T>(-0.2033*me_c);
    T pz =  static_cast<T>(94.5843*me_c);
    T ex =  static_cast<T>(39.8996*eref);
    T ey =  static_cast<T>(-29.2501*eref);
    T ez =  static_cast<T>(58.7720*eref);
    T bx =  static_cast<T>(44.3417*bref);
    T by =  static_cast<T>(15.5024*bref);
    T bz =  static_cast<T>(29.4024*bref);

    T exp = static_cast<T>(4.69766211952e-73*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL((exp-res)/exp, tolerance<T>());
}


//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_5 )
{
    breit_wheeler_engine_prod_5<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_5 )
{
    breit_wheeler_engine_prod_5<float>();
}

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


//Test pair production rates (generic)
template <typename T>
void breit_wheeler_engine_prod_7()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(61019.1*me_c);
    T py =  static_cast<T>(-24359.3*me_c);
    T pz =  static_cast<T>(65116.2*me_c);
    T ex =  static_cast<T>(69942.0*eref);
    T ey =  static_cast<T>(38024.7*eref);
    T ez =  static_cast<T>(-43604.1*eref);
    T bx =  static_cast<T>(-26990.0*bref);
    T by =  static_cast<T>(58267.8*bref);
    T bz =  static_cast<T>(-63485.8*bref);

    T exp = static_cast<T>(7.63488202211*rateref);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz})*light_speed, chi_photon(px,py,pz,ex,ey,ez,bx,by,bz) );

    BOOST_CHECK_SMALL((exp-res)/exp, tolerance<T>());
}

//Test pair production rates (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_7 )
{
    breit_wheeler_engine_prod_7<double>();
}

//Test pair production rates (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_single_7 )
{
    breit_wheeler_engine_prod_7<float>();
}

//Test evolve_opt_depth_and_determine_event (generic)
template <typename T>
void breit_wheeler_engine_detopt_1()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(149.825*me_c);
    T py =  static_cast<T>(933.115*me_c);
    T pz =  static_cast<T>(-538.195*me_c);
    T ex =  static_cast<T>(931.686*eref);
    T ey =  static_cast<T>(-861.074*eref);
    T ez =  static_cast<T>(944.652*eref);
    T bx =  static_cast<T>(531.406*bref);
    T by =  static_cast<T>(670.933*bref);
    T bz =  static_cast<T>(660.057*bref);

    T initial_optical_depth = static_cast<T>(1.0);
    T optical_depth = initial_optical_depth;
    T dt = static_cast<T>(0.01*dtref);

    T exp_rate = static_cast<T>(1.50648551484*rateref);

    bool has_event_happend;
    T dt_prod;

    bw_engine.compute_dN_dt_lookup_table(nullptr);

    std::tie(has_event_happend, dt_prod) =
        bw_engine.evolve_opt_depth_and_determine_event
        (px, py, pz, ex, ey, ez, bx, by, bz, dt, optical_depth);

    T exp_opt_depth = initial_optical_depth - dt*exp_rate;

    BOOST_CHECK_EQUAL(has_event_happend, false);
    BOOST_CHECK_EQUAL(dt_prod, static_cast<T>(0.0));
    BOOST_CHECK_SMALL((optical_depth-exp_opt_depth)/exp_opt_depth,
        tolerance<T>());
}


//Test evolve_opt_depth_and_determine_event (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_detopt_double_1 )
{
    breit_wheeler_engine_detopt_1<double>();
}

//Test evolve_opt_depth_and_determine_event (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_detopt_single_1 )
{
    breit_wheeler_engine_detopt_1<float>();
}


//Test evolve_opt_depth_and_determine_event (generic)
template <typename T>
void breit_wheeler_engine_detopt_2()
{
    auto bw_engine = get_bw_stl_set_lambda<T>(390109317);

    T px =  static_cast<T>(149.825*me_c);
    T py =  static_cast<T>(933.115*me_c);
    T pz =  static_cast<T>(-538.195*me_c);
    T ex =  static_cast<T>(931.686*eref);
    T ey =  static_cast<T>(-861.074*eref);
    T ez =  static_cast<T>(944.652*eref);
    T bx =  static_cast<T>(531.406*bref);
    T by =  static_cast<T>(670.933*bref);
    T bz =  static_cast<T>(660.057*bref);

    T initial_optical_depth = static_cast<T>(1.0e-3);
    T optical_depth = initial_optical_depth;
    T dt = static_cast<T>(0.01*dtref);

    T exp_rate = static_cast<T>(1.50648551484*rateref);

    bool has_event_happend;
    T dt_prod;

    bw_engine.compute_dN_dt_lookup_table(nullptr);

    std::tie(has_event_happend, dt_prod) =
        bw_engine.evolve_opt_depth_and_determine_event
        (px, py, pz, ex, ey, ez, bx, by, bz, dt, optical_depth);

    T exp_opt_depth = initial_optical_depth - dt*exp_rate;
    T exp_dt_prod = initial_optical_depth/exp_rate ;

    BOOST_CHECK_EQUAL(has_event_happend, true);
    BOOST_CHECK_SMALL(exp_dt_prod, dt_prod);
    BOOST_CHECK_SMALL((optical_depth-exp_opt_depth)/exp_opt_depth,
        tolerance<T>());
}


//Test evolve_opt_depth_and_determine_event (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_detopt_double_2 )
{
    breit_wheeler_engine_detopt_2<double>();
}

//Test evolve_opt_depth_and_determine_event (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_detopt_single_2 )
{
    breit_wheeler_engine_detopt_2<float>();
}

//Test evolve_opt_depth_and_determine_event (generic)
template <typename T>
void breit_wheeler_engine_detopt_3()
{
    breit_wheeler_engine_ctrl<T> bw_ctrl;
    bw_ctrl.chi_phot_min =  static_cast<T>(0.0001);

    auto bw_engine = get_bw_stl_set_lambda<T>(390109317, bw_ctrl);

    T px =  static_cast<T>(149.825*me_c);
    T py =  static_cast<T>(933.115*me_c);
    T pz =  static_cast<T>(-538.195*me_c);
    T ex =  static_cast<T>(931.686*eref);
    T ey =  static_cast<T>(-861.074*eref);
    T ez =  static_cast<T>(944.652*eref);
    T bx =  static_cast<T>(531.406*bref);
    T by =  static_cast<T>(670.933*bref);
    T bz =  static_cast<T>(660.057*bref);

    T initial_optical_depth = static_cast<T>(1.0e-3);
    T optical_depth = initial_optical_depth;
    T dt = static_cast<T>(0.01*dtref);

    T exp_rate = static_cast<T>(1.50648551484*rateref);

    bool has_event_happend;
    T dt_prod;

    bw_engine.compute_dN_dt_lookup_table(nullptr);

    std::tie(has_event_happend, dt_prod) =
    bw_engine.evolve_opt_depth_and_determine_event
    (px, py, pz, ex, ey, ez, bx, by, bz, dt, optical_depth);

    T exp_opt_depth = initial_optical_depth - dt*exp_rate;
    T exp_dt_prod = initial_optical_depth/exp_rate ;

    BOOST_CHECK_EQUAL(has_event_happend, true);
    BOOST_CHECK_SMALL(exp_dt_prod, dt_prod);
    BOOST_CHECK_SMALL((optical_depth-exp_opt_depth)/exp_opt_depth,
    tolerance<T>());
}


//Test evolve_opt_depth_and_determine_event (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_detopt_double_3 )
{
    breit_wheeler_engine_detopt_3<double>();
}

//Test evolve_opt_depth_and_determine_event (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_detopt_single_3 )
{
    breit_wheeler_engine_detopt_3<float>();
}
