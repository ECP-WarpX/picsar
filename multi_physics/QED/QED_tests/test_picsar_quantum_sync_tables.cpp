//####### Test module for Quantum Synchrotron tables ###########################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/tables"

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp>

#include <vector>
#include <algorithm>
#include <array>

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-3;
const double double_small = 1e-20;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-2;
const float float_small = 1e-10;


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

const double chi_min = 0.001;
const double chi_max = 1000;
const int how_many = 500;
const int how_many_frac = 500;
const int frac_first = 300;
const double frac_min = 1e-12;
const double frac_switch = 5e-2;

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
auto get_em_table()
{
    const auto params =
        photon_emission_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max),
            static_cast<RealType>(frac_min),
            how_many, how_many_frac};

    return photon_emission_lookup_table<RealType, VectorType>{params};
}

template <typename RealType, typename VectorType>
auto get_tailopt_em_table()
{
    const auto params =
        tailopt_photon_emission_lookup_table_params<RealType>{
            static_cast<RealType>(chi_min),
            static_cast<RealType>(chi_max),
            static_cast<RealType>(frac_min),
            static_cast<RealType>(frac_switch),
            how_many, how_many_frac, frac_first};

    return tailopt_photon_emission_lookup_table<RealType, VectorType>{params};
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

     for (int i = 0 ; i < static_cast<int>(coords.size()); ++i){
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

    const RealType ye_ext_o0 = alpha*chi_min;
    const RealType ye_ext_o1 = alpha*chi_max;
    const RealType ye0 = alpha*x0;
    const RealType ye1 = alpha*x1;
    const RealType ye2 = alpha*x2;

    const auto xxs = std::array<RealType, 5>
        {xo0, x0, x1, x2, xo1};

    const auto exp_ext = std::array<RealType, 5>
        {ye_ext_o0, ye0, ye1, ye2, ye_ext_o1};

    const auto is_out = std::array<bool, 5>
        {true, false, false, false, true};

    for(int i = 0 ; i < static_cast<int>(xxs.size()) ; ++i){
        const RealType res = table.interp(xxs[i]);
        bool flag_out = false;
        const RealType res2 = table.interp(xxs[i], &flag_out);
        BOOST_CHECK_EQUAL(flag_out, is_out[i]);
        BOOST_CHECK_EQUAL(res, res2);

        const RealType expect = exp_ext[i];

        if(i != 0)
            BOOST_CHECK_SMALL((res-expect)/expect, tolerance<RealType>());
        else
            BOOST_CHECK_SMALL((res-expect), tolerance<RealType>());
    }

    const auto table_view = table.get_view();

    for(int i = 0 ; i < static_cast<int>(xxs.size()) ; ++i){
        BOOST_CHECK_EQUAL(table_view.interp(xxs[i]), table.interp(xxs[i]));
    }
}

// ***Test Quantum Synchrotron dndt table
BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_dndt_table)
{
    check_dndt_table<double, std::vector<double>>();
    check_dndt_table<float, std::vector<float>>();
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

// ***Test Quantum Synchrotron dndt table serialization
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_dndt_table_serialization)
{
    check_dndt_table_serialization<float, std::vector<float>>();
}

template <typename RealType, typename VectorType>
void check_photon_emission_table()
{
    auto table = get_em_table<RealType, VectorType>();
    BOOST_CHECK_EQUAL(table.is_init(),false);

    auto coords = table.get_all_coordinates();

    BOOST_CHECK_EQUAL(coords.size(),how_many*how_many_frac);

    const RealType log_chi_min = log(chi_min);
    const RealType log_chi_max = log(chi_max);
    const RealType log_frac_min = log(frac_min);

     for (int i = 0 ; i < how_many*how_many_frac; ++i){
         auto res_1 = coords[i][0];
         auto res_2 = coords[i][1];
         const auto ii = i/how_many_frac;
         const auto jj = i%how_many_frac;
         auto expected_1 = static_cast<RealType>(
             exp(log_chi_min +ii*(log_chi_max-log_chi_min)/(how_many-1)));
        auto expected_2 = static_cast<RealType>(
             expected_1*exp(log_frac_min +jj*(0.0-log_frac_min)/(how_many_frac-1)));

         BOOST_CHECK_SMALL((res_1-expected_1)/expected_1, tolerance<RealType>());
         if(expected_2 != static_cast<RealType>(0.0))
            BOOST_CHECK_SMALL((res_2-expected_2)/expected_2, tolerance<RealType>());
        else
            BOOST_CHECK_SMALL((res_2-expected_2), tolerance<RealType>());
     }

    auto vals = VectorType(coords.size());

    auto functor = [=](std::array<RealType,2> x){
        return static_cast<RealType>(1.0*pow(x[1]/x[0], 4.0 + log10(x[0])));};

    auto inverse_functor = [=](std::array<RealType,2> x){
            return static_cast<RealType>(1.0*pow((1-x[1]), 1.0/(4.0 + log10(x[0]))));};

    std::transform(coords.begin(), coords.end(), vals.begin(),functor);

    bool result = table.set_all_vals(vals);
    BOOST_CHECK_EQUAL(result,true);
    BOOST_CHECK_EQUAL(table.is_init(),true);


    auto table_2 = get_em_table<RealType, VectorType>();
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
            auto expected = inverse_functor(std::array<RealType,2>{rxx, rr})*xx;
            if(expected < small<RealType>())
                BOOST_CHECK_SMALL((res-expected)/expected, tolerance<RealType>());
            else
                BOOST_CHECK_SMALL((res-expected), tolerance<RealType>());

        }
    }

    const auto table_view = table.get_view();

    for(auto xx : xxs){
        for (auto r : rrs)
            BOOST_CHECK_EQUAL(table_view.interp(xx,r ), table.interp(xx, r));
    }
}

// ***Test Quantum Synchrotron photon emission table
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_photon_emission_table)
{
    check_photon_emission_table<double, std::vector<double>>();
    check_photon_emission_table<float, std::vector<float>>();
}


template <typename RealType, typename VectorType>
void check_photon_emission_table_serialization()
{
    const auto params =
        photon_emission_lookup_table_params<RealType>{
            static_cast<RealType>(0.1),static_cast<RealType>(10.0),
            static_cast<RealType>(1e-5), 3, 3};

    auto table = photon_emission_lookup_table<
        RealType, VectorType>{params, {1.,2.,3.,4.,5.,6.,7.,8.,9.}};

    auto raw_data = table.serialize();
    auto new_table =
        photon_emission_lookup_table<RealType, VectorType>{raw_data};

    BOOST_CHECK_EQUAL(new_table.is_init(), true);
    BOOST_CHECK_EQUAL(new_table == table, true);
}


// ***Test Quantum Synchrotron tail-optimized photon emission table serialization
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_photon_emission_table_serialization)
{
    check_photon_emission_table_serialization<double, std::vector<double>>();
    check_photon_emission_table_serialization<float, std::vector<float>>();
}

template <typename RealType, typename VectorType>
void check_tailopt_photon_emission_table()
{
    auto table = get_tailopt_em_table<RealType, VectorType>();
    BOOST_CHECK_EQUAL(table.is_init(),false);

    auto coords = table.get_all_coordinates();

    BOOST_CHECK_EQUAL(coords.size(),how_many*how_many_frac);

    const auto lin_functor = detail::LinFunctor<RealType>{
        how_many,
        static_cast<RealType>(std::log(chi_min)),
        static_cast<RealType>(std::log(chi_max))};

    const auto tailopt_functor = detail::TailOptFunctor<RealType>{
        how_many_frac,
        frac_first,
        static_cast<RealType>(std::log(frac_min)),
        static_cast<RealType>(std::log(1.0)),
        static_cast<RealType>(std::log(frac_switch))};

    for (int i = 0 ; i < how_many*how_many_frac; ++i){

        auto res_1 = coords[i][0];
        auto res_2 = coords[i][1];
        const auto ii = i/how_many_frac;
        const auto jj = i%how_many_frac;

        auto expected_1 = std::exp(lin_functor(ii));
        auto expected_2 = std::exp(tailopt_functor(jj))*expected_1;

        BOOST_CHECK_SMALL((res_1-expected_1)/expected_1, tolerance<RealType>());
        if(expected_2 != static_cast<RealType>(0.0))
            BOOST_CHECK_SMALL((res_2-expected_2)/expected_2, tolerance<RealType>());
        else
            BOOST_CHECK_SMALL((res_2-expected_2), tolerance<RealType>());
    }

    auto vals = VectorType(coords.size());

    auto functor = [=](std::array<RealType,2> x){
        return static_cast<RealType>(1.0*pow(x[1]/x[0], 4.0 + log10(x[0])));};

    auto inverse_functor = [=](std::array<RealType,2> x){
            return static_cast<RealType>(1.0*pow((1-x[1]), 1.0/(4.0 + log10(x[0]))));};

    std::transform(coords.begin(), coords.end(), vals.begin(),functor);

    bool result = table.set_all_vals(vals);
    BOOST_CHECK_EQUAL(result,true);
    BOOST_CHECK_EQUAL(table.is_init(),true);

    auto table_2 = get_tailopt_em_table<RealType, VectorType>();
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
            {0.0, 1.0e-4, 1.0e-3, 1.0e-2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};

    for (const auto xx : xxs){
        for (const auto rr : rrs){
            auto res = table.interp(xx, rr);
            auto rxx = xx;
            if(rxx < chi_min) rxx = chi_min;
            if(rxx > chi_max) rxx = chi_max;
            auto expected = inverse_functor(std::array<RealType,2>{rxx, rr})*xx;
            if(expected < small<RealType>())
                BOOST_CHECK_SMALL((res-expected)/expected, tolerance<RealType>());
            else
                BOOST_CHECK_SMALL((res-expected), tolerance<RealType>());

        }
    }

    const auto table_view = table.get_view();

    for(auto xx : xxs){
        for (auto r : rrs)
            BOOST_CHECK_EQUAL(table_view.interp(xx,r ), table.interp(xx, r));
    }

}

// ***Test Quantum Synchrotron tail-optimized photon emission table
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_tailopt_photon_emission_table)
{
    check_tailopt_photon_emission_table<double, std::vector<double>>();
    check_tailopt_photon_emission_table<float, std::vector<float>>();
}


template <typename RealType, typename VectorType>
void check_tailopt_photon_emission_table_serialization()
{
    const auto params =
        tailopt_photon_emission_lookup_table_params<RealType>{
            static_cast<RealType>(0.1),static_cast<RealType>(10.0),
            static_cast<RealType>(1e-5), static_cast<RealType>(1e-1),
            3, 3, 2};

    auto table = tailopt_photon_emission_lookup_table<
        RealType, VectorType>{params, {1.,2.,3.,4.,5.,6.,7.,8.,9.}};

    auto raw_data = table.serialize();
    auto new_table =
        tailopt_photon_emission_lookup_table<RealType, VectorType>{raw_data};

    BOOST_CHECK_EQUAL(new_table.is_init(), true);
    BOOST_CHECK_EQUAL(new_table == table, true);
}


// ***Test Quantum Synchrotron tail-optimized photon emission table serialization
BOOST_AUTO_TEST_CASE( picsar_quantum_sync_tailopt_photon_emission_table_serialization)
{
    check_tailopt_photon_emission_table_serialization<double, std::vector<double>>();
    check_tailopt_photon_emission_table_serialization<float, std::vector<float>>();
}
