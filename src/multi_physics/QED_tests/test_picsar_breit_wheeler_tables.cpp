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
const double double_tolerance = 5.0e-5;
const double double_small = 1e-30;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-3;
const float float_small = 1e-20;


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

const double chi_min = 0.001;
const double chi_max = 1000;
const int how_many = 256;
const int how_many_frac = 256;

template <typename RealType, typename VectorType>
auto get_table()
{
    const auto params =
        dndt_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max), how_many};

    return dndt_lookup_table<RealType, VectorType>{params};
}

template <typename RealType, typename VectorType>
auto get_pair_table()
{
    const auto params =
        pair_prod_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max),
            how_many, how_many_frac};

    return pair_prod_lookup_table_logchi_linfrac<RealType, VectorType>{params};
}


// ------------- Tests --------------

template <typename RealType, typename VectorType>
void check_dndt_table()
{
    auto table = get_table<RealType, VectorType>();
    BOOST_CHECK_EQUAL(table.is_init(),false);

    VectorType coords = table.get_all_coordinates();

    BOOST_CHECK_EQUAL(coords.size(),how_many);

    const RealType log_chi_min = log(chi_min);
    const RealType log_chi_max = log(chi_max);

     for (int i = 0 ; i < coords.size(); ++i){
         auto res = coords[i];
         auto expected = static_cast<RealType>(
             exp(log_chi_min + i*(log_chi_max-log_chi_min)/(how_many-1)));
         BOOST_CHECK_SMALL((res-expected)/expected, tolerance<RealType>());
     }

    auto vals = VectorType{coords};

    const RealType alpha = 3.0;

    std::transform(coords.begin(), coords.end(), vals.begin(),
        [=](RealType x){return alpha*x;});

    bool result = table.set_all_vals(vals);
    BOOST_CHECK_EQUAL(result,true);
    BOOST_CHECK_EQUAL(table.is_init(),true);

    auto table_2 = get_table<RealType, VectorType>();
    BOOST_CHECK_EQUAL(table_2 == table, false);
    BOOST_CHECK_EQUAL(table == table, true);
    BOOST_CHECK_EQUAL(table_2 == table_2, true);

    const RealType xo0 = chi_min*0.1;
    const RealType xo1 = chi_max*10;
    const RealType x0 = chi_min;
    const RealType x1 = (chi_max+chi_min)*0.5642 + chi_min;
    const RealType x2 = chi_max;

    const RealType ye_app_o0 = dndt_approx_left<RealType>(xo0);
    const RealType ye_app_o1 = dndt_approx_right<RealType>(xo1);
    const RealType ye0 = alpha*x0;
    const RealType ye1 = alpha*x1;
    const RealType ye2 = alpha*x2;

    const auto xxs = std::array<RealType, 5>
        {xo0, x0, x1, x2, xo1};

    const auto exp_app = std::array<RealType, 5>
        {ye_app_o0, ye0, ye1, ye2, ye_app_o1};

    const auto is_out = std::array<bool, 5>
        {true, false, false, false, true};

    for(int i = 0 ; i < xxs.size() ; ++i){
        const RealType res = table.interp(xxs[i]);
        bool flag_out = false;
        const RealType res2 = table.interp_flag_out(xxs[i], flag_out);
        BOOST_CHECK_EQUAL(flag_out, is_out[i]);
        BOOST_CHECK_EQUAL(res, res2);

        const RealType expect = exp_app[i];

        if(i != 0)
            BOOST_CHECK_SMALL((res-expect)/expect, tolerance<RealType>());
        else
            BOOST_CHECK_SMALL((res-expect), tolerance<RealType>());
    }

    const auto table_view = table.get_view();

    for(int i = 0 ; i < xxs.size() ; ++i){
        BOOST_CHECK_EQUAL(table_view.interp(xxs[i]), table.interp(xxs[i]));
    }
}

// ***Test Breit Wheeler dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_table)
{
    check_dndt_table<double, std::vector<double>>();
    check_dndt_table<float, std::vector<float>>();
}

template <typename RealType>
void check_dndt_table_out_approx()
{
    const auto left_chi = std::array<RealType,4>{0.01,0.05,0.1,0.2};
    const auto left_sol = std::array<RealType,4>{
        2.9069621438923337e-117,
        1.2969667695320204e-24,
        4.944416339700773e-13,
        3.0528686912504136e-07};

    for(int i = 0; i < 4; i++){
        const auto left_res = dndt_approx_left(left_chi[i]);
        if(left_sol[i] < small<RealType>()){
            BOOST_CHECK_SMALL(
                left_res,small<RealType>());
        }
        else{
            BOOST_CHECK_SMALL(
                (left_res - left_sol[i])/left_sol[i],tolerance<RealType>());
        }
    }

}

// ***Test Breit Wheeler dndt table out range approximation
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_table_out_approx)
{
    check_dndt_table_out_approx<double>();
    check_dndt_table_out_approx<float>();
}

template <typename RealType, typename VectorType>
void check_dndt_table_serialization()
{
    const auto params =
        dndt_lookup_table_params<RealType>{
            static_cast<RealType>(0.1),
            static_cast<RealType>(10.0), 3};

    auto table = dndt_lookup_table<
        RealType, VectorType>{params, {1.,2.,3.}};

    auto raw_data = table.serialize();
    auto new_table = dndt_lookup_table<
        RealType, VectorType>{raw_data};

    BOOST_CHECK_EQUAL(new_table.is_init(), true);
    BOOST_CHECK_EQUAL(new_table == table, true);
}

// ***Test Breit Wheeler dndt table serialization
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_table_serialization)
{
    check_dndt_table_serialization<double, std::vector<double>>();
    check_dndt_table_serialization<float, std::vector<float>>();
}

template <typename RealType, typename VectorType>
void check_pair_production_table()
{
    auto table = get_pair_table<RealType, VectorType>();
    BOOST_CHECK_EQUAL(table.is_init(),false);

    auto coords = table.get_all_coordinates();

    BOOST_CHECK_EQUAL(coords.size(),how_many*how_many_frac);

    const RealType log_chi_min = log(chi_min);
    const RealType log_chi_max = log(chi_max);

     for (int i = 0 ; i < how_many*how_many_frac; ++i){
         auto res_1 = coords[i][0];
         auto res_2 = coords[i][1];
         const auto ii = i/how_many_frac;
         const auto jj = i%how_many_frac;
         auto expected_1 = static_cast<RealType>(
             exp(log_chi_min +ii*(log_chi_max-log_chi_min)/(how_many-1)));
        auto expected_2 = static_cast<RealType>(
             0.0 +jj*0.5/(how_many_frac-1));

         BOOST_CHECK_SMALL((res_1-expected_1)/expected_1, tolerance<RealType>());
         if(expected_2 != static_cast<RealType>(0.0))
            BOOST_CHECK_SMALL((res_2-expected_2)/expected_2, tolerance<RealType>());
        else
            BOOST_CHECK_SMALL((res_2-expected_2), tolerance<RealType>());
     }

    auto vals = VectorType(coords.size());

    auto functor = [=](std::array<RealType,2> x){
        return static_cast<RealType>(0.5*pow(2*x[1], 8.0 + log(x[0])));};

    auto inverse_functor = [=](std::array<RealType,2> x){
            return static_cast<RealType>(0.5*pow(2*x[1], 1.0/(8.0 + log(x[0]))));};

    std::transform(coords.begin(), coords.end(), vals.begin(),functor);

    bool result = table.set_all_vals(vals);
    BOOST_CHECK_EQUAL(result,true);
    BOOST_CHECK_EQUAL(table.is_init(),true);


    auto table_2 = get_pair_table<RealType, VectorType>();
    BOOST_CHECK_EQUAL(table_2 == table, false);
    BOOST_CHECK_EQUAL(table == table, true);
    BOOST_CHECK_EQUAL(table_2 == table_2, true);


    const RealType xo0 = chi_min*0.1;
    const RealType xo1 = chi_max*10;
    const RealType x0 = chi_min;
    const RealType x1 = (chi_max+chi_min)*0.5642 + chi_min;
    const RealType x2 = chi_max;

    const auto xxs = std::array<RealType, 5>
        {xo0, x0, x1, x2, xo1};

    const auto rrs = std::array<RealType, 11>
            {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};

    for (const auto xx : xxs){
        for (const auto rr : rrs){
            auto res = table.interp(xx, rr);
            auto rxx = xx;
            if(rxx < chi_min) rxx = chi_min;
            if(rxx > chi_max) rxx = chi_max;
            auto expected = inverse_functor(std::array<RealType,2>{rxx, rr/2})*xx;
            if(rr >= 0.5) expected = xx - expected;

            if(expected != static_cast<RealType>(0.0))
                BOOST_CHECK_SMALL((res-expected)/expected, tolerance<RealType>());
            else
                BOOST_CHECK_SMALL((res-expected), tolerance<RealType>());

        }
    }
    const auto table_view = table.get_view();

    const auto ff = std::array<double,4>{0.0, 0.1, 0.5, 0.99};

    for(int i = 0 ; i < xxs.size() ; ++i){
        for (auto f : ff)
            BOOST_CHECK_EQUAL(table_view.interp(xxs[i],f ), table.interp(xxs[i], f));
    }
}

// ***Test Breit Wheeler dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_pair_production_table)
{
    check_pair_production_table<double, std::vector<double>>();
    check_pair_production_table<float, std::vector<float>>();
    check_pair_production_table<double, std::vector<double>>();
    check_pair_production_table<float, std::vector<float>>();
}


template <typename RealType, typename VectorType>
void check_pair_production_table_serialization()
{
    const auto params =
        pair_prod_lookup_table_params<RealType>{
            static_cast<RealType>(0.1),static_cast<RealType>(10.0),
            3, 3};

    auto table = pair_prod_lookup_table_logchi_linfrac<
        RealType, VectorType>{params, {1.,2.,3.,4.,5.,6.,7.,8.,9.}};

    auto raw_data = table.serialize();
    auto new_table =
        pair_prod_lookup_table_logchi_linfrac<RealType, VectorType>{raw_data};

    BOOST_CHECK_EQUAL(new_table.is_init(), true);
    BOOST_CHECK_EQUAL(new_table == table, true);
}


// ***Test Breit Wheeler pair production table serialization
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_pair_production_table_serialization)
{
    check_pair_production_table_serialization<double, std::vector<double>>();
    check_pair_production_table_serialization<float, std::vector<float>>();
}
