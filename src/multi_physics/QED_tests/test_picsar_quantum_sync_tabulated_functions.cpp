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
const double double_tolerance = 3.0e-2;
//Tolerance for single precision calculations
const float float_tolerance = 3.0e-2;

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

// ------------- Tests --------------

template <typename RealType>
void check_int_K_5_3_replacement()
{
    const auto cases = std::array<std::pair<RealType,RealType>,5>{
        std::make_pair<RealType,RealType>( 1e-5, 4629.2044114881355),
        std::make_pair<RealType,RealType>( 1e-4, 995.9088308508012),
        std::make_pair<RealType,RealType>( 1e-2, 44.497250411420913),
        std::make_pair<RealType,RealType>( 1, 0.651422815355309),
        std::make_pair<RealType,RealType>( 10, 1.9223826430338323e-05)};

    for (const auto cc : cases)
    {
        const auto res = inner_integral(cc.first);
        const auto sol = cc.second;
        BOOST_CHECK_SMALL((res-sol)/sol,  tolerance<RealType>());
    }

}

// ***Test replacement of the integral of K_5_3
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_int_K_5_3_replacement)
{
    check_int_K_5_3_replacement<double>();
    check_int_K_5_3_replacement<float>();
}


template <typename RealType>
void check_dndt_table()
{
    const auto cases = std::array<std::pair<double,double>,9>{
        std::make_pair( 0.0001, 0.0002164863577428999),
        std::make_pair( 0.001, 0.0021630710451102635),
        std::make_pair( 0.01, 0.02145778613250966),
        std::make_pair( 0.1, 0.20141650057288696),
        std::make_pair( 1.0 , 1.5508709239783094),
        std::make_pair( 10.0, 9.170292626506058),
        std::make_pair( 100.0, 46.02341774244706),
        std::make_pair( 1000.0, 217.8438638968691),
        std::make_pair( 10000.0, 1015.6105987224346)};

    for (const auto cc : cases)
    {
        const auto res = compute_G_function(static_cast<RealType>(cc.first));
            BOOST_CHECK_SMALL((res - static_cast<RealType>(cc.second))/
                static_cast<RealType>(cc.second), tolerance<RealType>());
    }
}

// ***Test Quantum Synchrotron dndt table
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_dndt_G_function)
{
    check_dndt_table<double>();
    check_dndt_table<float>();
}
