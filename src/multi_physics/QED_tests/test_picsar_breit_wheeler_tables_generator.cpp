//####### Test module for Breit-Wheeler tables generator #######################

//Define Module name
 #define BOOST_TEST_MODULE "phys/breit_wheeler/tables_generator"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <vector>
#include <algorithm>
#include <array>

#include "breit_wheeler_engine_tables.hpp"
#include "breit_wheeler_engine_tables_generator.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-2;
const double double_small = 1e-16;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-2;
const float float_small = 1e-8;


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

template <typename T>
T constexpr small()
{
    if(std::is_same<T,float>::value)
        return float_small;
    else
        return double_small;
}

const double chi_min = 0.1;
const double chi_max = 100;
const int how_many = 47;
const int how_many_frac = 29;

// ------------- Tests --------------

// ***Test Breit Wheeler dndt table generation

template <typename RealType, typename VectorType>
void check_dndt_table_generation()
{
    const auto params =
        dndt_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max), how_many};

    auto table = dndt_lookup_table<RealType, VectorType>{params};

    table.generate();

    const auto chi_T_vector = std::vector<std::array<RealType,2>>{
            std::array<RealType,2>{0.001, 0.0},
            std::array<RealType,2>{0.01, 3.535432626057024e-117},
            std::array<RealType,2>{0.1, 5.92605753015639e-13},
            std::array<RealType,2>{1.0, 0.014135754351952334},
            std::array<RealType,2>{10.0, 0.10848609251601983},
            std::array<RealType,2>{100.0, 0.07489672111818155},
            std::array<RealType,2>{1000.0, 0.03728739404639084}};

    for (const auto chi_T : chi_T_vector){
        bool is_out = false;
        const auto res = table.interp(chi_T[0], &is_out);
        const auto exp = chi_T[1];

        if(exp > small<RealType>()){
            BOOST_CHECK_SMALL((res-exp)/exp,tolerance<RealType>());
        }
        else{
            BOOST_CHECK_SMALL(res,small<RealType>());
        }
    }

}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_table_generation)
{
    check_dndt_table_generation<double, std::vector<double>>();
    check_dndt_table_generation<float, std::vector<float>>();
}

// *******************************

// ***Test Breit Wheeler pair production table generation

template <typename RealType, typename VectorType>
void check_pair_prod_table_generation()
{
    const auto params =
        pair_prod_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max), how_many,how_many_frac};

    auto table = pair_prod_lookup_table<RealType, VectorType>{params};

    table.generate();



    BOOST_CHECK_EQUAL(true,true);
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_pair_prod_table_generation)
{
    check_pair_prod_table_generation<double, std::vector<double>>();
    check_pair_prod_table_generation<float, std::vector<float>>();
}

// *******************************
