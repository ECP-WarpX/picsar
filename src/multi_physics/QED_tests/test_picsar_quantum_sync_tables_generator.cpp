//####### Test module for  Quantum Synchrotron tables generator ################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/tables"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <vector>
#include <algorithm>
#include <array>

#include "quantum_sync_engine_tables.hpp"
#include "quantum_sync_engine_tables_generator.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-2;
const double double_small = 1e-16;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-2;
const float float_small = 1e-8;


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

const double chi_min = 0.01;
const double chi_max = 100;
const int how_many = 29;
const int how_many_frac = 47;
const double frac_min = 1e-12;

// ------------- Tests --------------

// ***Test Quantum Synchrotron dndt table generation

template <typename RealType, typename VectorType>
void check_dndt_table_generation()
{
    const auto params =
        dndt_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max), how_many};

    auto table = dndt_lookup_table<RealType, VectorType>{params};

    table.generate();

    const auto chi_G_vector = std::vector<std::array<RealType,2>>{
            std::array<RealType,2>{0.001, 0.02145778613250966}, // out of table: 0.01 is used
            std::array<RealType,2>{0.01, 0.02145778613250966},
            std::array<RealType,2>{0.1, 0.20141650057288696},
            std::array<RealType,2>{1.0, 1.5508709239783094},
            std::array<RealType,2>{10.0, 9.170292626506058},
            std::array<RealType,2>{100.0, 46.02341774244706},
            std::array<RealType,2>{1000.0,46.02341774244706}}; // out of table: 100 is used

    for (const auto chi_G : chi_G_vector){
        bool is_out = false;
        const auto res = table.interp(chi_G[0]);
        const auto exp = chi_G[1];

        BOOST_CHECK_EQUAL(is_out, (chi_G[0] < chi_min ) || (chi_G[0] > chi_max) );

        if(exp > small<RealType>()){
            BOOST_CHECK_SMALL((res-exp)/exp,tolerance<RealType>());
        }
        else{
            BOOST_CHECK_SMALL(res,small<RealType>());
        }
    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_dndt_table_generation)
{
    check_dndt_table_generation<double, std::vector<double>>();
    check_dndt_table_generation<float, std::vector<float>>();
}

// *******************************

// ***Test Quantum Synchrotron photon emission table generation

template <typename RealType, typename VectorType>
void check_photon_emission_table_generation()
{
    const auto params =
        photon_emission_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max),
            static_cast<RealType>(frac_min),
            how_many,how_many_frac};

    auto table = photon_emission_lookup_table<RealType, VectorType>{params};

    table.generate();

    const auto chi_chi_P_vector = std::vector<std::array<RealType,3>>{
        std::array<RealType,3>{ 1.0 , 0.3 , 0.1345442547767904 },
        std::array<RealType,3>{ 1.0 , 0.5 , 0.5 },
        std::array<RealType,3>{ 1.0 , 0.7 , 0.8654557452351763 },
        std::array<RealType,3>{ 10.0 , 1.0 , 0.07575595226082421 },
        std::array<RealType,3>{ 10.0 , 3.0 , 0.3096952419129383 },
        std::array<RealType,3>{ 10.0 , 5.0 , 0.5 },
        std::array<RealType,3>{ 10.0 , 7.0 , 0.6903047580870976 },
        std::array<RealType,3>{ 10.0 , 9.0 , 0.924244047739035 },
        std::array<RealType,3>{ 100.0 , 10.0 , 0.16909092328399775 },
        std::array<RealType,3>{ 100.0 , 30.0 , 0.36646213315061366 },
        std::array<RealType,3>{ 100.0 , 50.0 , 0.5000000000000001 },
        std::array<RealType,3>{ 100.0 , 70.0 , 0.633537866849343 },
        std::array<RealType,3>{ 100.0 , 90.0 , 0.8309090767162309 },
        // OUT OF TABLE, USE VALUES FOR CHI = 100
        std::array<RealType,3>{ 1000.0 , 500.0 , 0.5000000000000001 },
        std::array<RealType,3>{ 1000.0 , 700.0 , 0.633537866849343 },
        std::array<RealType,3>{ 1000.0 , 900.0 , 0.8309090767162309 }
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

    BOOST_CHECK_EQUAL(true,true);
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_photon_emission_table_generation)
{
    check_photon_emission_table_generation<double, std::vector<double>>();
    check_photon_emission_table_generation<float, std::vector<float>>();
}

// *******************************
