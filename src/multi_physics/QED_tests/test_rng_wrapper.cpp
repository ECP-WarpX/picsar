//####### Test module for rng_wrapper ####################################

//Define Module name
 #define BOOST_TEST_MODULE "RNG wrapper"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

#include <cstdint>
#include <random>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS
#include "rng_wrapper.hpp"

//Only include Kokkos if it is supported
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT
    #include <Kokkos_Core.hpp>
    #include <Kokkos_Random.hpp>
    #include <Kokkos_DualView.hpp>
#endif

using namespace picsar::multi_physics;



// ------------- Tests --------------

//***STL***

//***Constructors

//Test STL rng_wrapper constructors(double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_constructors_double )
{
    int64_t seed = 2391892344079;
    std::mt19937_64 rng{2391892344079};

    stl_rng_wrapper wrp1{seed};
    stl_rng_wrapper wrp2(move(rng));

    BOOST_CHECK_EQUAL( wrp1.unf<double>(0.0,1.0), wrp2.unf<double>(0.0,1.0));
}

//Test STL rng_wrapper constructors(single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_constructors_single )
{
    int64_t seed = 2391892344079;
    std::mt19937_64 rng{2391892344079};

    stl_rng_wrapper wrp1{seed};
    stl_rng_wrapper wrp2(move(rng));

    BOOST_CHECK_EQUAL( wrp1.unf<float>(0.0f,1.0f), wrp2.unf<float>(0.0f,1.0f));
}

//***Uniform distribution

//Test STL rng_wrapper unf (double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_unf_double )
{
    int64_t seed = 2391892344079;
    stl_rng_wrapper wrp{seed};
    size_t how_many = 10000;
    double a = -7.0;
    double b = 11.1;

    for (size_t d = 0; d <= how_many; ++d){
        double qq = wrp.unf<double> (a,b);
        BOOST_TEST( qq >= a);
        BOOST_TEST( qq < b);
    }
}

//Test STL rng_wrapper unf (single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_unf_single )
{
    int64_t seed = 2391892344079;
    stl_rng_wrapper wrp{seed};
    size_t how_many = 10000;
    float a = -7.0f;
    float b = 11.1f;

    for (size_t d = 0; d <= how_many; ++d){
        float qq = wrp.unf<float>(a,b);
        BOOST_TEST( qq >= a);
        BOOST_TEST( qq < b);
    }
}

//Test STL rng_wrapper exp (double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_exp_double )
{
    int64_t seed = 2391892344079;
    stl_rng_wrapper wrp{seed};
    size_t how_many = 10000;
    double l = 1.0f;

    for (size_t d = 0; d <= how_many; ++d){
        double qq = wrp.exp<double>(l);
        BOOST_TEST( qq >= 0.0);
    }
}

//Test STL rng_wrapper exp (single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_exp_single )
{
    int64_t seed = 2391892344079;
    stl_rng_wrapper wrp{seed};
    size_t how_many = 10000;
    float l = 1.0f;

    for (size_t d = 0; d <= how_many; ++d){
        double qq = wrp.exp<float>(l);
        BOOST_TEST( qq >= 0.0f);
    }
}

//***Kokkos***
//Build ONLY if Kokkos is enabled
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

//Test Kokkos rng_wrapper unf (double precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_unf_double_1 )
{
    //Initialize Kokkos passing argc and argv provided through Boost test suite
    Kokkos::initialize(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv);

    //Creates a pointer to a pool of generators
    auto pool = Kokkos::Random_XorShift1024_Pool<>{89139811};

    //We want to fill a 1D array with numbers extracted from a uniform
    //distribution between a and b
    size_t size = 100000;
    double a = -7.0;
    double b = 11.0;

    //A view is the way Kokkos manages multi-dimensional arrays
    //A dual view manages mirroring between a view living in the device memory
    //and a view living in the host memory.
    Kokkos::DualView<double*> vals("A bunch of random numbers", size);

    //I create my wrapper
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrp{pool};

    //I create a lambda, which should be passed to Kokkos::parallel_for
    auto ffunc =
    KOKKOS_LAMBDA (int i){
        vals.d_view(i) = wrp.unf(a,b);
    };

    Kokkos::parallel_for(size, ffunc);
    Kokkos::fence();

    //Copy from device to host
    Kokkos::deep_copy(vals.h_view,vals.d_view);

    Kokkos::finalize();

    for(size_t i = 0; i < size; ++i){
        BOOST_TEST( vals.h_view(i) >= a);
        BOOST_TEST( vals.h_view(i) < b);
    }

}

//Test Kokkos rng_wrapper exp (double precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_exp_double_1 )
{
    //Initialize Kokkos passing argc and argv provided through Boost test suite
    Kokkos::initialize(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv);

    //Creates a pointer to a pool of generators
    auto pool = Kokkos::Random_XorShift1024_Pool<>{89139811};

    //We want to fill a 1D array with numbers extracted from an exp distribution
    size_t size = 100000;

    //A view is the way Kokkos manages multi-dimensional arrays
    //A dual view manages mirroring between a view living in the device memory
    //and a view living in the host memory.
    Kokkos::DualView<double*> vals("A bunch of random numbers", size);

    //I create my wrapper
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrp{pool};

    //I create a lambda, which should be passed to Kokkos::parallel_for
    auto ffunc =
    KOKKOS_LAMBDA (int i){
        vals.d_view(i) = wrp.exp(1.0);
    };

    Kokkos::parallel_for(size, ffunc);
    Kokkos::fence();

    //Copy from device to host
    Kokkos::deep_copy(vals.h_view,vals.d_view);

    Kokkos::finalize();

    for(size_t i = 0; i < size; ++i){
        BOOST_TEST( vals.h_view(i) >= 0.0);
    }

}


//Test Kokkos rng_wrapper unf (single precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_unf_single_1 )
{
    //Initialize Kokkos passing argc and argv provided through Boost test suite
    Kokkos::initialize(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv);

    //Creates a pointer to a pool of generators
    auto pool = Kokkos::Random_XorShift1024_Pool<>{89139811};

    //We want to fill a 1D array with numbers extracted from a uniform
    //distribution between a and b
    size_t size = 100000;
    float a = -7.0;
    float b = 11.0;

    //A view is the way Kokkos manages multi-dimensional arrays
    //A dual view manages mirroring between a view living in the device memory
    //and a view living in the host memory.
    Kokkos::DualView<float*> vals("A bunch of random numbers", size);

    //I create my wrapper
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrp{pool};

    //I create a lambda, which should be passed to Kokkos::parallel_for
    auto ffunc =
    KOKKOS_LAMBDA (int i){
        vals.d_view(i) = wrp.unf<float>(a,b);
    };

    Kokkos::parallel_for(size, ffunc);
    Kokkos::fence();

    //Copy from device to host
    Kokkos::deep_copy(vals.h_view,vals.d_view);

    Kokkos::finalize();

    for(size_t i = 0; i < size; ++i){
        BOOST_TEST( vals.h_view(i) >= a);
        BOOST_TEST( vals.h_view(i) < b);
    }

}

//Test Kokkos rng_wrapper exp (single precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_exp_single_1 )
{
    //Initialize Kokkos passing argc and argv provided through Boost test suite
    Kokkos::initialize(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv);

    //Creates a pointer to a pool of generators
    auto pool = Kokkos::Random_XorShift1024_Pool<>{89139811};

    //We want to fill a 1D array with numbers extracted from an exp distribution
    size_t size = 100000;

    //A view is the way Kokkos manages multi-dimensional arrays
    //A dual view manages mirroring between a view living in the device memory
    //and a view living in the host memory.
    Kokkos::DualView<double*> vals("A bunch of random numbers", size);

    //I create my wrapper
    kokkos_rng_wrapper<Kokkos::Random_XorShift1024_Pool<>> wrp{pool};

    //I create a lambda, which should be passed to Kokkos::parallel_for
    auto ffunc =
    KOKKOS_LAMBDA (int i){
        vals.d_view(i) = wrp.exp(1.0f);
    };

    Kokkos::parallel_for(size, ffunc);
    Kokkos::fence();

    //Copy from device to host
    Kokkos::deep_copy(vals.h_view,vals.d_view);

    Kokkos::finalize();

    for(size_t i = 0; i < size; ++i){
        BOOST_TEST( vals.h_view(i) >= 0.0f);
    }

}

#endif

//***Exponential distribution
