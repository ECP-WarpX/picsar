//####### Test module for  Quantum Synchrotron tables generator ################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/tables"

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <vector>
#include <algorithm>
#include <array>

#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp>
#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp>

//Tolerance for double precision calculations
const double double_tolerance = 3.0e-2;
const double double_small = 1e-12;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-2;
const float float_small = 1e-6;


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
const int how_many = 73;
const int how_many_frac = 75;
const double frac_min = 1e-6;

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
        const auto res = table.interp(chi_G[0],&is_out);
        const auto exp = chi_G[1];

        BOOST_CHECK_EQUAL(is_out, (chi_G[0] < RealType(chi_min) ) || (chi_G[0] > RealType(chi_max)) );

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
        const auto res = table.interp(chi_chi_P[0], RealType(1.) - chi_chi_P[2], &is_out);
        const auto exp = chi_chi_P[1];

        BOOST_CHECK_EQUAL(is_out, (chi_chi_P[0] < RealType(chi_min) ) || (chi_chi_P[0] > RealType(chi_max)) );

        if(exp > small<RealType>()){
            BOOST_CHECK_SMALL((res-exp)/exp,tolerance<RealType>());
        }
        else{
            BOOST_CHECK_SMALL(res,small<RealType>());
        }
    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_photon_emission_table_generation)
{
    check_photon_emission_table_generation<double, std::vector<double>>();
    check_photon_emission_table_generation<float, std::vector<float>>();
}

// *******************************
