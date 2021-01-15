//####### Test module for quantum sync tabulated functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/tabulated_functions"

 #include<array>
 #include<utility>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp>

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-7;
const double double_small = 1e-22;

//Tolerance for single precision calculations
const float float_tolerance = 5.0e-4;
const float float_small = 1e-16;

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

// ***Test replacement the integral of K(5/3, x)

template <typename RealType>
void check_int_K_5_3_replacement()
{
    const auto cases = std::array<std::pair<RealType,RealType>,5>{
        std::make_pair<RealType,RealType>( 1e-5, 4629.2044114881355),
        std::make_pair<RealType,RealType>( 1e-4, 995.9088308508012),
        std::make_pair<RealType,RealType>( 1e-2, 44.497250411420913),
        std::make_pair<RealType,RealType>( 1, 0.651422815355309),
        std::make_pair<RealType,RealType>( 10, 1.9223826430338323e-05)};

    for (const auto& cc : cases)
    {
        const auto res = inner_integral(cc.first);
        const auto sol = cc.second;
        BOOST_CHECK_SMALL((res-sol)/sol,  tolerance<RealType>());
    }

}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_int_K_5_3_replacement)
{
    check_int_K_5_3_replacement<double>();
    check_int_K_5_3_replacement<float>();
}

// *******************************

// ***Test G function

template <typename RealType>
void check_G_function()
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

    for (const auto& cc : cases)
    {
        const auto res = compute_G_function(static_cast<RealType>(cc.first));
            BOOST_CHECK_SMALL((res - static_cast<RealType>(cc.second))/
                static_cast<RealType>(cc.second), tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_G_function)
{
    check_G_function<double>();
    check_G_function<float>();
}

// *******************************

// ***Test cumulative probability distribution

template <typename RealType>
void check_cumulative_prob_1()
{
    const auto chiphot_chipart_res = std::array<std::array<double,3>,50>{
        std::array<double,3>{ 0.0 , 0.0 , 0.0},
        std::array<double,3>{ 0.001, 0.0 , 0.0},
        std::array<double,3>{ 0.001, 1e-06, 0.8571415045564768},
        std::array<double,3>{ 0.001, 0.0003 , 1.0},
        std::array<double,3>{ 0.001, 0.0005, 1.0},
        std::array<double,3>{ 0.001, 0.0007, 1.0},
        std::array<double,3>{ 0.001, 0.000999, 1.0},
        std::array<double,3>{ 0.001, 0.001, 1.0 },
        std::array<double,3>{ 0.01, 0.0, 0.0},
        std::array<double,3>{ 0.01, 1e-05, 0.480547190762555},
        std::array<double,3>{ 0.01, 0.003, 1.0},
        std::array<double,3>{ 0.01, 0.005, 1.0},
        std::array<double,3>{ 0.01, 0.007, 1.0},
        std::array<double,3>{ 0.01, 0.00999, 1.0},
        std::array<double,3>{ 0.01, 0.01, 1.0},
        std::array<double,3>{ 0.1, 0.0,  0.0},
        std::array<double,3>{ 0.1, 0.0001, 0.246635968071763},
        std::array<double,3>{ 0.1, 0.03, 0.9956252388535388},
        std::array<double,3>{ 0.1, 0.05, 0.9999555773722274},
        std::array<double,3>{ 0.1, 0.06999999999999999, 0.9999999974091136},
        std::array<double,3>{ 0.1, 0.0999, 1.0},
        std::array<double,3>{ 0.1, 0.1, 1.0},
        std::array<double,3>{ 1.0, 0.0, 0.0},
        std::array<double,3>{ 1.0, 0.001, 0.1498504179995893 },
        std::array<double,3>{ 1.0, 0.3, 0.8646332229371583 },
        std::array<double,3>{ 1.0, 0.5, 0.9509588179865093 },
        std::array<double,3>{ 1.0, 0.7, 0.9904326213067777 },
        std::array<double,3>{ 1.0, 0.999, 1.0 },
        std::array<double,3>{ 1.0, 1 , 1.0 },
        std::array<double,3>{ 10.0, 0.0, 0.0 },
        std::array<double,3>{ 10.0, 0.01, 0.11782827804283932 },
        std::array<double,3>{ 10.0, 3.0, 0.7376845730958391 },
        std::array<double,3>{ 10.0, 5.0, 0.8452405292352005 },
        std::array<double,3>{ 10.0, 7.0, 0.9209032377215558 },
        std::array<double,3>{ 10.0, 9.99, 0.9999999999999206 },
        std::array<double,3>{ 10.0, 10.0, 1.0 },
        std::array<double,3>{ 100.0, 0.0, 0.0 },
        std::array<double,3>{ 100.0, 0.1, 0.10901296836655985 },
        std::array<double,3>{ 100.0, 30.0, 0.6941324168235529 },
        std::array<double,3>{ 100.0, 50.0, 0.8017131199967799 },
        std::array<double,3>{ 100.0, 70.0, 0.8807461242142075 },
        std::array<double,3>{ 100.0, 99.9, 0.9999995295418699 },
        std::array<double,3>{ 100.0, 100.0, 1.0 },
        std::array<double,3>{ 1000.0, 0.0 ,  0.0 },
        std::array<double,3>{ 1000.0, 1.0 ,  0.10690852121367304 },
        std::array<double,3>{ 1000.0, 300.0, 0.6831879431840134 },
        std::array<double,3>{ 1000.0, 500.0, 0.7903662468832419 },
        std::array<double,3>{ 1000.0, 700.0, 0.8696031647277049 },
        std::array<double,3>{ 1000.0, 999.0, 0.9993578963869831 },
        std::array<double,3>{ 1000.0, 1000.0, 1.0 }};

    for (const auto cc : chiphot_chipart_res)
    {
        const auto chi_part = static_cast<RealType>(cc[0]);
        const auto chi_phot = static_cast<RealType>(cc[1]);
        const auto expected = static_cast<RealType>(cc[2]);
        const auto res = compute_cumulative_prob(
            chi_part, std::vector<RealType>{chi_phot})[0];
        if(expected < small<RealType>()){
            BOOST_CHECK_SMALL( (res - expected), small<RealType>());
        }else{
            BOOST_CHECK_SMALL((res - expected)/expected, tolerance<RealType>());
        }

    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_cumulative_prob_1)
{
    check_cumulative_prob_1<double>();
    check_cumulative_prob_1<float>();
}

// *******************************

// ***Test cumulative probability distribution (vector of inputs)

template <typename RealType>
void check_cumulative_prob_vec()
{
    const RealType chipart = 10.0;
    const auto chiphot = std::vector<RealType>{
        0.0, 1e-3, 1e-2, 1e-1, 1, 2, 3, 4, 5,
        6, 7, 8, 9, 9.9, 9.99, 9.999, 10.0};
    const auto expected = std::vector<RealType>{0.0,
        0.0547191031257403, 0.11782827804283932,
        0.25304594229286587, 0.5329925008830805,
        0.65721264691597, 0.7376845730958391,
        0.7974768847498939, 0.8452405292352005,
        0.8854544941535698,0.9209032377215558,
        0.9534724891171927,0.9837225233708096,
        0.9999975010317529,0.9999999999999206,
        1.0000000000002132,1.0};
    const auto res = compute_cumulative_prob(chipart, chiphot);

    for (int i = 0 ; i < static_cast<int>(expected.size()); ++i)
    {
        if(expected[i] < small<RealType>()){
            BOOST_CHECK_SMALL( (res[i] - expected[i]), small<RealType>());
        }else{
            BOOST_CHECK_SMALL((res[i] - expected[i])/expected[i], tolerance<RealType>());
        }

    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_cumulative_prob_vec)
{
    check_cumulative_prob_vec<double>();
    check_cumulative_prob_vec<float>();
}

// *******************************
