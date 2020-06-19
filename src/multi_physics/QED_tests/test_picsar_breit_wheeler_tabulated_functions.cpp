//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/breit_wheeler/tabulated_functions"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

 #include<array>
 #include<utility>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "breit_wheeler_engine_tabulated_functions.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-8;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-3;

using namespace picsar::multi_physics::phys::breit_wheeler;

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
void check_dndt_table()
{
    const auto cases = std::array<std::pair<double,double>,20>{
        std::make_pair( 0.0001, 0.0),
        std::make_pair( 0.00026366508987303583, 0.0),
        std::make_pair( 0.0006951927961775605, 0.0),
        std::make_pair( 0.0018329807108324356, 0.0),
        std::make_pair( 0.004832930238571752, 5.369549410542618e-241),
        std::make_pair( 0.012742749857031334, 2.9894334826923954e-92),
        std::make_pair( 0.03359818286283781, 7.742890804075623e-36),
        std::make_pair( 0.08858667904100823, 1.911813697083578e-14),
        std::make_pair( 0.23357214690901212, 2.43856256298753e-06),
        std::make_pair( 0.615848211066026, 0.002783057482755572),
        std::make_pair( 1.623776739188721, 0.0374292083884315),
        std::make_pair( 4.281332398719396, 0.08995768126397806),
        std::make_pair( 11.288378916846883, 0.10884131132341256),
        std::make_pair( 29.763514416313132, 0.09912946656214293),
        std::make_pair( 78.47599703514607, 0.07986454178500298),
        std::make_pair( 206.913808111479, 0.06088073823190407),
        std::make_pair( 545.5594781168514, 0.04521754955598924),
        std::make_pair( 1438.44988828766, 0.033161084780181496),
        std::make_pair( 3792.690190732246, 0.02416544874852762),
        std::make_pair( 10000.0, 0.017553195930973545)};

    for (const auto cc : cases)
    {
        const auto res = compute_T_function(static_cast<RealType>(cc.first));
        if(cc.second < tolerance<RealType>()){
            BOOST_CHECK_SMALL( res, tolerance<RealType>());
        }else{
            BOOST_CHECK_SMALL((res - static_cast<RealType>(cc.second))/
                static_cast<RealType>(cc.second), tolerance<RealType>());
        }

    }
}

// ***Test Breit Wheeler dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_T_function)
{
    check_dndt_table<double>();
    check_dndt_table<float>();
}