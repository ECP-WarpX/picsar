//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/breit_wheeler/tables"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <vector>
#include <algorithm>
#include <array>

#include "breit_wheeler_engine_tables.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-5;

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

template <typename RealType, typename VectorType, dndt_table_type TableType>
auto get_table()
{
    const double chi_min = 0.001;
    const double chi_max = 1000;
    const int how_many = 100;

    const auto params =
        dndt_lookup_table_params<RealType>{chi_min, chi_max, how_many};

    return dndt_lookup_table<
        RealType, VectorType, dndt_table_type::logchi>{params};

}
// ------------- Tests --------------

template <typename RealType, typename VectorType, dndt_table_type TableType>
void check_dndt_table()
{
    auto table = get_table<RealType, VectorType, TableType>();
    BOOST_CHECK_EQUAL(0,0);
}

// ***Test Breit Wheeler dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_table)
{
    check_dndt_table<double, std::vector<double>, dndt_table_type::logchi>();
    check_dndt_table<float, std::vector<float>, dndt_table_type::logchi>();
}
