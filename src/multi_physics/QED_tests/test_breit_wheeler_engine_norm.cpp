//####### Test module for Breit Wheeler engine (SI units) ######################

//Define Module name
 #define BOOST_TEST_MODULE "Nonliner Breit Wheeler engine (normalized units)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include <utility>

 //Include Boost unit tests library & library for floating point comparison
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>

//SI units are used for this test
#define PXRMP_USE_NORMALIZED_UNITS
#include "breit_wheeler_engine.hpp"

#include "chi_functions.hpp"
#include "vec_functions.hpp"

using namespace picsar::multi_physics;

//________________________________

//Helper function
template<typename REAL, typename WHATEVER>
breit_wheeler_engine<REAL, stl_rng_wrapper> get_bw_stl_set_lambda(uint64_t seed, WHATEVER lambda,
breit_wheeler_engine_ctrl<REAL> bw_ctrl = breit_wheeler_engine_ctrl<REAL>())
{
    stl_rng_wrapper wrap{seed};
    auto bw_engine =  breit_wheeler_engine<REAL, stl_rng_wrapper>{std::move(wrap), 1.0, bw_ctrl};
    bw_engine.set_lambda(static_cast<REAL>(lambda));
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

//Special tolerance for a specific case (check carefully)

//Special Tolerance for double precision calculations
const double double_spec_tolerance = 2.0e-2;

//Special Tolerance for single precision calculations
const float float_spec_tolerance = 2.0e-2;

//Special Templated tolerance
template <typename T>
T spec_tolerance()
{
    if(std::is_same<T,float>::value)
        return float_spec_tolerance;
    else
        return double_spec_tolerance;
}

//Test get/set lambda for breit_wheeler_engine generic
template <typename T>
void breit_wheeler_engine_gs()
{
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));
    BOOST_CHECK_EQUAL( static_cast<T>(800.0*si_nanometer), bw_engine.get_lambda());
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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));
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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(0.0);
    T py =  static_cast<T>(0.0);
    T pz =  static_cast<T>(0.0);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(-861.074);
    T ez =  static_cast<T>(944.652);
    T bx =  static_cast<T>(531.406);
    T by =  static_cast<T>(670.933);
    T bz =  static_cast<T>(660.057);
    T lambda = static_cast<T>(800. * si_nanometer);

    //T exp = static_cast<T>(0.0);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(149.825);
    T py =  static_cast<T>(0.0);
    T pz =  static_cast<T>(0.0);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(0.0);
    T ez =  static_cast<T>(0.0);
    T bx =  static_cast<T>(0.0);
    T by =  static_cast<T>(0.0);
    T bz =  static_cast<T>(0.0);
    T lambda = static_cast<T>(800. * si_nanometer);

    //T exp = static_cast<T>(0.0);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(6.81696);
    T py =  static_cast<T>(9.68933);
    T pz =  static_cast<T>(2.81229);
    T ex =  static_cast<T>(-4.89986);
    T ey =  static_cast<T>(-9.65535);
    T ez =  static_cast<T>(3.69471);
    T bx =  static_cast<T>(8.89549);
    T by =  static_cast<T>(-5.46574);
    T bz =  static_cast<T>(-6.75393);
    T lambda = static_cast<T>(800. * si_nanometer);

    //T exp = static_cast<T>(0.0);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(149.825);
    T py =  static_cast<T>(933.115);
    T pz =  static_cast<T>(-538.195);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(-861.074);
    T ez =  static_cast<T>(944.652);
    T bx =  static_cast<T>(531.406);
    T by =  static_cast<T>(670.933);
    T bz =  static_cast<T>(660.057);
    T lambda = static_cast<T>(800. * si_nanometer);

    T exp = static_cast<T>(1.50648551484);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

    std::cerr << chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) <<" " << exp << " " << res << std::endl;

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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(-44.4546);
    T py =  static_cast<T>(-0.2033);
    T pz =  static_cast<T>(94.5843);
    T ex =  static_cast<T>(39.8996);
    T ey =  static_cast<T>(-29.2501);
    T ez =  static_cast<T>(58.7720);
    T bx =  static_cast<T>(44.3417);
    T by =  static_cast<T>(15.5024);
    T bz =  static_cast<T>(29.4024);
    T lambda = static_cast<T>(800. * si_nanometer);

    //T exp = static_cast<T>(4.69766211952e-73);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

    BOOST_CHECK_SMALL(res, tolerance<T>());
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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(-965.61);
    T py =  static_cast<T>(-3975.11);
    T pz =  static_cast<T>(6917.22);
    T ex =  static_cast<T>(11.17);
    T ey =  static_cast<T>(-2117.72);
    T ez =  static_cast<T>(-1407.19);
    T bx =  static_cast<T>( 6259.79);
    T by =  static_cast<T>(7557.54);
    T bz =  static_cast<T>(773.11);
    T lambda = static_cast<T>(800. * si_nanometer);

    T exp = static_cast<T>(3.51855878777);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));
    T px =  static_cast<T>(61019.1);
    T py =  static_cast<T>(-24359.3);
    T pz =  static_cast<T>(65116.2);
    T ex =  static_cast<T>(69942.0);
    T ey =  static_cast<T>(38024.7);
    T ez =  static_cast<T>(-43604.1);
    T bx =  static_cast<T>(-26990.0);
    T by =  static_cast<T>(58267.8);
    T bz =  static_cast<T>(-63485.8);
    T lambda = static_cast<T>(800. * si_nanometer);

    T exp = static_cast<T>(7.63488202211);
    T res = bw_engine.compute_dN_dt(norm(vec3<T>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

    BOOST_CHECK_SMALL((exp-res)/exp, spec_tolerance<T>());
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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(149.825);
    T py =  static_cast<T>(933.115);
    T pz =  static_cast<T>(-538.195);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(-861.074);
    T ez =  static_cast<T>(944.652);
    T bx =  static_cast<T>(531.406);
    T by =  static_cast<T>(670.933);
    T bz =  static_cast<T>(660.057);

    T initial_optical_depth = static_cast<T>(1.0);
    T optical_depth = initial_optical_depth;
    T dt = static_cast<T>(0.01);

    T exp_rate = static_cast<T>(1.50648551484);

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
    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer));

    T px =  static_cast<T>(149.825);
    T py =  static_cast<T>(933.115);
    T pz =  static_cast<T>(-538.195);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(-861.074);
    T ez =  static_cast<T>(944.652);
    T bx =  static_cast<T>(531.406);
    T by =  static_cast<T>(670.933);
    T bz =  static_cast<T>(660.057);

    T initial_optical_depth = static_cast<T>(1.0e-3);
    T optical_depth = initial_optical_depth;
    T dt = static_cast<T>(0.01);

    T exp_rate = static_cast<T>(1.50648551484);

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

    auto bw_engine = get_bw_stl_set_lambda<T>
        (390109317, static_cast<T>(800.0*si_nanometer), bw_ctrl);

    T px =  static_cast<T>(149.825);
    T py =  static_cast<T>(933.115);
    T pz =  static_cast<T>(-538.195);
    T ex =  static_cast<T>(931.686);
    T ey =  static_cast<T>(-861.074);
    T ez =  static_cast<T>(944.652);
    T bx =  static_cast<T>(531.406);
    T by =  static_cast<T>(670.933);
    T bz =  static_cast<T>(660.057);

    T initial_optical_depth = static_cast<T>(1.0e-3);
    T optical_depth = initial_optical_depth;
    T dt = static_cast<T>(0.01);

    T exp_rate = static_cast<T>(1.50648551484);

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
