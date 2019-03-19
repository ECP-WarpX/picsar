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

//Test STL rng_wrapper constructors generic
template<typename T>
void  rng_stl_wrapper_constructors(int64_t seed)
{
    std::mt19937_64 rng{seed};

    stl_rng_wrapper wrp1{seed};
    stl_rng_wrapper wrp2(move(rng));

    BOOST_CHECK_EQUAL( wrp1.unf<T>(0.0,1.0), wrp2.unf<T>(0.0,1.0));
}

//Test STL rng_wrapper constructors(double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_constructors_double )
{
    rng_stl_wrapper_constructors<double>(2391892344079);
}

//Test STL rng_wrapper constructors(single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_constructors_single )
{
    rng_stl_wrapper_constructors<float>(2391892344079);
}

//***Uniform distribution

//Test STL rng_wrapper unf generic
template<typename T>
void  rng_stl_wrapper_unf(int64_t seed)
{
    stl_rng_wrapper wrp{seed};
    size_t how_many = 10000;
    T a = static_cast<T>(-7.0);
    T b = static_cast<T>(11.1);

    for (size_t d = 0; d <= how_many; ++d){
        T qq = wrp.unf<T> (a,b);
        BOOST_TEST( qq >= a);
        BOOST_TEST( qq < b);
    }
}

//Test STL rng_wrapper unf (double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_unf_double )
{
    rng_stl_wrapper_unf<double>(2391892344079);
}

//Test STL rng_wrapper unf (single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_unf_single )
{
    rng_stl_wrapper_unf<float>(2391892344079);
}


//Test STL rng_wrapper exp generic
template<typename T>
void  rng_stl_wrapper_exp(int64_t seed)
{
    stl_rng_wrapper wrp{seed};
    size_t how_many = 10000;
    T l = static_cast<T>(1.0);

    for (size_t d = 0; d <= how_many; ++d){
        T qq = wrp.exp<T>(l);
        BOOST_TEST( qq >= 0.0);
    }
}


//Test STL rng_wrapper exp (double precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_exp_double )
{
    rng_stl_wrapper_exp<double>(2391892344079);
}

//Test STL rng_wrapper exp (single precision)
BOOST_AUTO_TEST_CASE( rng_stl_wrapper_exp_single )
{
    rng_stl_wrapper_exp<float>(2391892344079);
}

//***Kokkos***
//Build ONLY if Kokkos is enabled
#ifdef PXRMP_BUILD_WITH_KOKKOS_SUPPORT

//Test Kokkos rng_wrapper generic
template<typename T>
void rng_kokkos_unf ()
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
    T a = -7.0;
    T b = 11.0;

    //A view is the way Kokkos manages multi-dimensional arrays
    //A dual view manages mirroring between a view living in the device memory
    //and a view living in the host memory.
    Kokkos::DualView<T*> vals("A bunch of random numbers", size);

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

    //Single test
    T val = wrp.unf(a,b);
    BOOST_TEST( val >= a);
    BOOST_TEST( val < b);

    Kokkos::finalize();

    for(size_t i = 0; i < size; ++i){
        BOOST_TEST( vals.h_view(i) >= a);
        BOOST_TEST( vals.h_view(i) < b);
    }

}

//Test Kokkos rng_wrapper unf (double precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_unf_double_1 )
{
    rng_kokkos_unf<double>();
}

//Test Kokkos rng_wrapper unf (single precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_unf_single_1 )
{
    rng_kokkos_unf<float>();
}


//Test Kokkos rng_wrapper exp generic
template<typename T>
void rng_kokkos_exp()
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
    Kokkos::DualView<T*> vals("A bunch of random numbers", size);

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

    //Single test
    BOOST_TEST( wrp.exp(static_cast<T>(1.0)) >= static_cast<T>(0.0));

    Kokkos::finalize();

    for(size_t i = 0; i < size; ++i){
        BOOST_TEST( vals.h_view(i) >= static_cast<T>(0.0));
    }

}


//Test Kokkos rng_wrapper exp (double precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_exp_double_1 )
{
    rng_kokkos_exp<double>();
}

//Test Kokkos rng_wrapper exp (single precision)
BOOST_AUTO_TEST_CASE( rng_kokkos_exp_single_1 )
{
    rng_kokkos_exp<float>();
}

#endif

//***Exponential distribution
