//####### Test module for breit wheeler tabulated functions ####################################

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
const double double_tolerance = 1.0e-7;
const double double_small = 1e-22;

//Tolerance for single precision calculations
const float float_tolerance = 5.0e-4;
const float float_small = 1e-16;

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

// ------------- Tests --------------

// ***Test T function

template <typename RealType>
void check_T_func()
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
        if(cc.second < small<RealType>()){
            BOOST_CHECK_SMALL( res, small<RealType>());
        }else{
            BOOST_CHECK_SMALL((res - static_cast<RealType>(cc.second))/
                static_cast<RealType>(cc.second), tolerance<RealType>());
        }

    }
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_T_function)
{
    check_T_func<double>();
    check_T_func<float>();
}

// *******************************

// ***Test cumulative probability distribution

template <typename RealType>
void check_cumulative_prob_1()
{
    const auto chiphot_chipart_res = std::array<std::array<double,3>,50>{
        std::array<double,3>{0.0, 0.0, 0.0},
        std::array<double,3>{0.001, 0.0, 0.0},
        std::array<double,3>{0.001, 1e-06, 0.0},
        std::array<double,3>{0.001, 0.0003, 0.0},
        std::array<double,3>{0.001, 0.0005, 0.5},
        std::array<double,3>{0.001, 0.0007, 1.0},
        std::array<double,3>{0.001, 0.000999, 1.0},
        std::array<double,3>{0.001, 0.001, 1.0},
        std::array<double,3>{0.01, 0.0, 0.0},
        std::array<double,3>{0.01, 1e-05, 0.0},
        std::array<double,3>{0.01, 0.003, 3.0227480641243037e-24},
        std::array<double,3>{0.01, 0.005, 0.5},
        std::array<double,3>{0.01, 0.006999999999999999, 1.0},
        std::array<double,3>{0.01, 0.00999, 0.9999999786081522},
        std::array<double,3>{0.01, 0.01, 1.0},
        std::array<double,3>{0.1, 0.0, 0.0},
        std::array<double,3>{0.1, 0.0001, 0.0},
        std::array<double,3>{0.1, 0.03, 0.0006311081101750777},
        std::array<double,3>{0.1, 0.05, 0.5000000000000001},
        std::array<double,3>{0.1, 0.06999999999999999, 0.9993688919143366},
        std::array<double,3>{0.1, 0.0999, 1.000000000000005},
        std::array<double,3>{0.1, 0.1, 1.0},
        std::array<double,3>{1, 0, 0.0},
        std::array<double,3>{1, 0.001, 1.4530789387331725e-293},
        std::array<double,3>{1, 0.3, 0.1345442547767904},
        std::array<double,3>{1, 0.5, 0.5},
        std::array<double,3>{1, 0.7, 0.8654557452351763},
        std::array<double,3>{1, 0.999, 1.0000000000006934},
        std::array<double,3>{1, 1.0, 1.0},
        std::array<double,3>{10.0, 0.0, 0.0},
        std::array<double,3>{10.0, 0.01, 3.981434709202223e-33},
        std::array<double,3>{10.0, 3.0, 0.3096952419129383},
        std::array<double,3>{10.0, 5.0, 0.5},
        std::array<double,3>{10.0, 7.0, 0.6903047580870976},
        std::array<double,3>{10.0, 9.99, 0.9999999999997944},
        std::array<double,3>{10.0, 10.0, 1.0},
        std::array<double,3>{100.0, 0.0, 0.0},
        std::array<double,3>{100.0, 0.1, 1.9014438671301059e-06},
        std::array<double,3>{100.0, 30.0, 0.36646213315061366},
        std::array<double,3>{100.0, 50.0, 0.5000000000000001},
        std::array<double,3>{100.0, 70.0, 0.633537866849343},
        std::array<double,3>{100.0, 99.9, 0.9999980985561322},
        std::array<double,3>{100.0, 100.0, 1.0},
        std::array<double,3>{1000.0, 0.0, 0.0},
        std::array<double,3>{1000.0, 1.0, 0.00249628362862492},
        std::array<double,3>{1000.0, 300.0, 0.37995481425306316},
        std::array<double,3>{1000.0, 500.0, 0.49999999999999994},
        std::array<double,3>{1000.0, 700.0, 0.6200451857470134},
        std::array<double,3>{1000.0, 999.0, 0.9975037163713809},
        std::array<double,3>{1000.0, 1000.0, 1.0}};

    for (const auto cc : chiphot_chipart_res)
    {
        const auto chi_phot = static_cast<RealType>(cc[0]);
        const auto chi_part = static_cast<RealType>(cc[1]);
        const auto expected = static_cast<RealType>(cc[2]);
        const auto res = compute_cumulative_prob(
            chi_phot, std::vector<RealType>{chi_part})[0];
        if(expected < small<RealType>()){
            BOOST_CHECK_SMALL( (res - expected), small<RealType>());
        }else{
            BOOST_CHECK_SMALL((res - expected)/expected, tolerance<RealType>());
        }

    }
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_cumulative_prob_1)
{
    check_cumulative_prob_1<double>();
    check_cumulative_prob_1<float>();
}

// *******************************

// ***Test cumulative probability distribution (vector of inputs)

template <typename RealType>
void check_cumulative_prob_vec()
{
    const RealType chiphot = 10.0;
    const auto chipart = std::vector<RealType>{
        0.0, 1e-3, 1e-2, 1e-1, 1, 2, 3, 4, 5,
        6, 7, 8, 9, 9.9, 9.99, 9.999, 10.0};
    const auto expected = std::vector<RealType>{
        0.0,3.45376351890413e-295,3.981434709202223e-33,1.2297633168391708e-05,
        0.07575595226082421,0.19941889404675717,0.3096952419129383,0.40793955045875263,
        0.5,0.592060449541101,0.6903047580870976,0.8005811059529608,0.924244047739035,
        0.9999877023668319,0.9999999999997944,0.9999999999999992,1.0};
    const auto res = compute_cumulative_prob(chiphot, chipart);

    for (int i = 0 ; i < expected.size(); ++i)
    {
        if(expected[i] < small<RealType>()){
            BOOST_CHECK_SMALL( (res[i] - expected[i]), small<RealType>());
        }else{
            BOOST_CHECK_SMALL((res[i] - expected[i])/expected[i], tolerance<RealType>());
        }

    }
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_cumulative_prob_vec)
{
    check_cumulative_prob_vec<double>();
    check_cumulative_prob_vec<float>();
}

// *******************************

