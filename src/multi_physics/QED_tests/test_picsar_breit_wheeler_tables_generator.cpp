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

        BOOST_CHECK_EQUAL(is_out, (chi_T[0] < chi_min ) || (chi_T[0] > chi_max) );

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

    const auto chi_chi_P_vector = std::vector<std::array<RealType,3>>{
        // OUT OF TABLE, USE VALUES FOR CHI = 0.1
        std::array<RealType,3>{ 0.01 , 0.00001 , 0.246635968071763 },
        std::array<RealType,3>{ 0.01 , 0.0001 , 0.5111462944662963 },
        std::array<RealType,3>{ 0.01 , 0.0010000000000000002 , 0.9029456666866598 },
        //_______________________________________
        std::array<RealType,3>{ 0.1 , 0.0001 , 0.246635968071763 },
        std::array<RealType,3>{ 0.1 , 0.001 , 0.5111462944662963 },
        std::array<RealType,3>{ 0.1 , 0.010000000000000002 , 0.9029456666866598 },
        std::array<RealType,3>{ 1.0 , 0.001 , 0.1498504179995893 },
        std::array<RealType,3>{ 1.0 , 0.01 , 0.31983031171402876 },
        std::array<RealType,3>{ 1.0 , 0.1 , 0.6537737809484806 },
        std::array<RealType,3>{ 1.0 , 0.2 , 0.7870086703462721 },
        std::array<RealType,3>{ 1.0 , 0.5 , 0.9509588179865093 },
        std::array<RealType,3>{ 10.0 , 0.01 , 0.11782827804283932 },
        std::array<RealType,3>{ 10.0 , 0.1 , 0.25304594229286587 },
        std::array<RealType,3>{ 10.0 , 1.0 , 0.5329925008830805 },
        std::array<RealType,3>{ 10.0 , 2.0 , 0.65721264691597 },
        std::array<RealType,3>{ 10.0 , 5.0 , 0.8452405292352005 },
        std::array<RealType,3>{ 10.0 , 9.0 , 0.9837225233708096 },
        std::array<RealType,3>{ 100.0 , 0.1 , 0.10901296836655985 },
        std::array<RealType,3>{ 100.0 , 1.0 , 0.2344249952917742 },
        std::array<RealType,3>{ 100.0 , 10.0 , 0.49689370771825486 },
        std::array<RealType,3>{ 100.0 , 20.0 , 0.6157367664025167 },
        std::array<RealType,3>{ 100.0 , 50.0 , 0.8017131199967799 },
        std::array<RealType,3>{ 100.0 , 90.0 , 0.9545809417165405 },
        // OUT OF TABLE, USE VALUES FOR CHI = 100
        std::array<RealType,3>{ 1000.0 , 1.0 , 0.10901296836655985 },
        std::array<RealType,3>{ 1000.0 , 10.0 , 0.2344249952917742 },
        std::array<RealType,3>{ 1000.0 , 100.0 , 0.49689370771825486 },
        std::array<RealType,3>{ 1000.0 , 200.0 , 0.6157367664025167 },
        std::array<RealType,3>{ 1000.0 , 500.0 , 0.8017131199967799 },
        std::array<RealType,3>{ 1000.0 , 900.0 , 0.9545809417165405 }
        //_______________________________________
        };

    for (const auto chi_chi_P : chi_chi_P_vector){
        bool is_out = false;
        const auto res = table.interp(chi_chi_P[0], chi_chi_P[2], &is_out);
        const auto exp = chi_chi_P[1];

        BOOST_CHECK_EQUAL(is_out, false);

        if(exp > small<RealType>()){
            BOOST_CHECK_SMALL((res-exp)/exp,tolerance<RealType>());
        }
        else{
            BOOST_CHECK_SMALL(res,small<RealType>());
        }
    }
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_pair_prod_table_generation)
{
    check_pair_prod_table_generation<double, std::vector<double>>();
    check_pair_prod_table_generation<float, std::vector<float>>();
}

// *******************************
