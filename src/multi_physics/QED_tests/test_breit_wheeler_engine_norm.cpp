//####### Test module for Breit Wheeler engine (SI units) ######################

//Define Module name
 #define BOOST_TEST_MODULE "Nonliner Breit Wheeler engine (SI)"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//SI units are used for this test
#define PXRMP_USE_NORMALIZED_UNITS
#include "breit_wheeler_engine.hpp"

#include "chi_functions.hpp"
#include "vec_functions.hpp"

//Only include Kokkos if it is supported
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT
    #include <Kokkos_Core.hpp>
    #include <Kokkos_Random.hpp>
    #include <Kokkos_DualView.hpp>
#endif

using namespace picsar::multi_physics;

//Helper function
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

//Kokkos version
template<typename REAL, typename WHATEVER>
breit_wheeler_engine
<REAL, kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>>>
get_bw_kokkos_set_lambda(int64_t seed, WHATEVER lambda)
{
    auto pool = Kokkos::Random_XorShift1024_Pool<>{seed};
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrap{pool};
    auto bw_engine =  breit_wheeler_engine
        <REAL, kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>>>{std::move(wrap)};
    bw_engine.set_lambda(static_cast<REAL>(lambda));
    return bw_engine;
}

#endif

//No Kokkos version
template<typename REAL, typename WHATEVER>
breit_wheeler_engine<REAL, stl_rng_wrapper> get_bw_stl_set_lambda(int64_t seed, WHATEVER lambda)
{
    stl_rng_wrapper wrap{seed};
    auto bw_engine =  breit_wheeler_engine<REAL, stl_rng_wrapper>{std::move(wrap)};
    bw_engine.set_lambda(static_cast<REAL>(lambda));
    return bw_engine;
}



// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-3;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-2;

//Test get/set lambda for breit_wheeler_engine (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_double_1 )
{
    auto bw_engine = get_bw_stl_set_lambda<double>
        (390109317, 800.0*si_nanometer);
    BOOST_CHECK_EQUAL( 800.0*si_nanometer, bw_engine.get_lambda());
}

//Test get/set lambda for breit_wheeler_engine (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_single_1 )
{
    auto bw_engine = get_bw_stl_set_lambda<float>
        (390109317, 800.0*si_nanometer);
    BOOST_CHECK_EQUAL( static_cast<float>(800.0*si_nanometer),
        bw_engine.get_lambda());
}

// ------------- optical depth --------------

//Test get new optical depth (double precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_opt_double_1 )
{
    auto bw_engine = get_bw_stl_set_lambda<double>
        (390109317, 800.0*si_nanometer);
    BOOST_TEST ( bw_engine.get_optical_depth() >= 0.0 );
}

//Test get new optical depth (single precision)
BOOST_AUTO_TEST_CASE( breit_wheeler_engine_opt_single_1 )
{
    auto bw_engine = get_bw_stl_set_lambda<float>
        (390109317, 800.0*si_nanometer);
    BOOST_TEST( bw_engine.get_optical_depth() >= 0.0 );
}

// ------------- Pair production rates --------------

//Test pair production rates (double precision)

BOOST_AUTO_TEST_CASE( breit_wheeler_engine_prod_double_1 )
{
    auto bw_engine = get_bw_stl_set_lambda<double>
        (390109317, 800.0*si_nanometer);

        double px = 149.825;
        double py = 933.115;
        double pz = -538.195;
        double ex = 931.686;
        double ey = -861.074;
        double ez = 944.652;
        double bx = 531.406;
        double by = 670.933;
        double bz = 660.057;
        double lambda = 800. * si_nanometer;

        double exp = 1.50648551484;

        double res = bw_engine.compute_dN_dt(norm(vec3<double>{px,py,pz}), chi_photon(px,py,pz,ex,ey,ez,bx,by,bz,lambda) );

    BOOST_TEST ( bw_engine.get_optical_depth() >= 0.0 );
}

// auto nn = [](double px,double py, double pz){return sqrt(px*px + py*py + pz*pz);};
// cout << "calc BW dN/dt (mom=[61019.1, -24359.3, 65116.2], EB=[69942.0, 38024.7, -43604.1, -26990.0, 58267.8, -63485.8], l = 800 nm, exp. 7.63488202211) : " << endl;
// cout << breit_wheeler_engine.get_total_pair_production_rate(nn(61019.1, -24359.3, 65116.2),
//     picsar::multi_physics::chi_photon_lambda({61019.1, -24359.3, 65116.2},{69942.0, 38024.7, -43604.1, -26990.0, 58267.8, -63485.8}, 0.8 * picsar::multi_physics::_um)) << endl;
// cout << "calc BW dN/dt (mom=[-965.61, -3975.11, 6917.22], EB=[11.17, -2117.72, -1407.19, 6259.79, 7557.54, 773.11], l = 800 nm, exp. 3.51855878777) : " << endl;
// cout << breit_wheeler_engine.get_total_pair_production_rate(nn(-965.61, -3975.11, 6917.22),
//     picsar::multi_physics::chi_photon_lambda({-965.61, -3975.11, 6917.22},{11.17, -2117.72, -1407.19, 6259.79, 7557.54, 773.11}, 0.8 * picsar::multi_physics::_um)) << endl;
// cout << "calc BW dN/dt (mom=[149.825, 933.115, -538.195], EB=[931.686, -861.074, 944.652, 531.406, 670.933, 660.057], l = 800 nm, exp. 1.50648551484) : " << endl;
// cout << breit_wheeler_engine.get_total_pair_production_rate(nn(149.825, 933.115, -538.195),
//     picsar::multi_physics::chi_photon_lambda({149.825, 933.115, -538.195},{931.686, -861.074, 944.652, 531.406, 670.933, 660.057}, 0.8 * picsar::multi_physics::_um)) << endl;
// cout << "calc BW dN/dt (mom=[-44.4546, -0.2033, 94.5843], EB=[39.8996, -29.2501, 58.7720, 44.3417, 15.5024, 29.4024], l = 800 nm, exp. 4.69766211952e-73) : " << endl;
// cout << breit_wheeler_engine.get_total_pair_production_rate(nn(-44.4546, -0.2033, 94.5843),
//     picsar::multi_physics::chi_photon_lambda({-44.4546, -0.2033, 94.5843},{39.8996, -29.2501, 58.7720, 44.3417, 15.5024, 29.4024}, 0.8 * picsar::multi_physics::_um)) << endl;
// cout << "calc BW dN/dt (mom=[6.81696,9.68933,2.81229], EB=[-4.89986,-9.65535,3.69471, 8.89549,-5.46574,-6.75393], l = 800 nm, exp. 0.0) : " << endl;
// cout << breit_wheeler_engine.get_total_pair_production_rate(nn(6.81696,9.68933,2.81229),
//     picsar::multi_physics::chi_photon_lambda({6.81696,9.68933,2.81229},{-4.89986,-9.65535,3.69471, 8.89549,-5.46574,-6.75393}, 0.8 * picsar::multi_physics::_um)) << endl;
