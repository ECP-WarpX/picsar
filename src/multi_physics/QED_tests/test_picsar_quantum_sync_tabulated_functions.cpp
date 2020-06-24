//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/tabulated_functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include<array>
 #include<utility>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "quantum_sync_engine_tabulated_functions.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-8;
const double double_small = 1e-30;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-3;
const float float_small = 1e-20;

using namespace picsar::multi_physics::phys::quantum_sync;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

template <typename T>
T constexpr small()
{
    if(std::is_same<T,float>::value)
        return float_small;
    else
        return double_small;
}

// ------------- Tests --------------

template <typename RealType>
void check_dndt_table()
{
    const auto cases = std::array<std::pair<double,double>,10>{
        std::make_pair( 1e-05, 4.074974691564334e-36),
        std::make_pair( 0.0001, 0.00021648693827504642),
        std::make_pair( 0.001, 0.002163071033312093),
        std::make_pair( 0.01, 0.021457689816994918),
        std::make_pair( 0.1, 0.20115778013982663),
        std::make_pair( 1.0 , 1.5053889224836552),
        std::make_pair( 10.0, 7.988917185870076),
        std::make_pair( 100.0, 36.82598230220699),
        std::make_pair( 1000.0, 152.5459946710829),
        std::make_pair( 10000.0, 518.3259462825962)};

    for (const auto cc : cases)
    {
        const auto res = compute_G_function(static_cast<RealType>(cc.first));
        if(cc.second < small<RealType>()){
            BOOST_CHECK_SMALL( res, small<RealType>());
        }else{
            BOOST_CHECK_SMALL((res - static_cast<RealType>(cc.second))/
                static_cast<RealType>(cc.second), tolerance<RealType>());
        }

    }
}

// ***Test Quantum Synchrotron dndt table
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_dndt_G_function)
{
    check_dndt_table<double>();
    check_dndt_table<float>();
}
