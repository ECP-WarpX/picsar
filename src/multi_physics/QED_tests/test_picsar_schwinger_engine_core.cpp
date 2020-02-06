//####### Test module for schwinger engine (core functions) ####################################

//Define Module name
 #define BOOST_TEST_MODULE "phys/schwinger"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <array>

#include "schwinger_pair_engine_core.hpp"

using namespace picsar::multi_physics::utils;

using namespace picsar::multi_physics::phys;

using namespace picsar::multi_physics::phys::schwinger;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-6;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-3;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//This wrapper is used for debug purposes to understand if the right
//functions are actually called
class wrapper_wrapper{
public:
    wrapper_wrapper(size_t seed):
        m_rng{seed}{};

    template <typename T>
    T unf(T a, T b)
    {
        auto val = m_rng.unf(a, b);
        m_flags.unf_calls++;
        m_flags.last_unf_call = val;
        return val;
    }

    template <typename T>
    T exp(T l)
    {
        auto val = m_rng.exp(l);
        m_flags.exp_calls++;
        m_flags.last_exp_call = val;
        return val;
    }

    template <typename T>
    size_t poisson(T l)
    {
        auto val = m_rng.poisson(l);
        m_flags.poisson_calls++;
        m_flags.last_poisson_call = val;
        return val;
    }

    template <typename T>
    T gaussian(T m, T s)
    {
        auto val = m_rng.gaussian(m, s);
        m_flags.gauss_calls++;
        m_flags.last_gauss_call = val;
        return val;
    }


    struct {
        int unf_calls = 0;
        double last_unf_call = 0.0;
        int exp_calls = 0;
        double last_exp_call = 0.0;
        int poisson_calls = 0;
        size_t last_poisson_call = 0;
        int gauss_calls = 0;
        double last_gauss_call = 0.0;
    } m_flags;

private:
    stl_rng_wrapper m_rng;
};


const double lambda = 800e-9;

// ------------- Tests --------------

std::vector<std::array<double,3>> E =
    {std::array<double,3>{3.57893317e+18,  5.46020302e+19, -9.93627604e+18}};

std::vector<std::array<double,3>> B =
    {std::array<double,3>{7.66560594e+10,  -1.58756483e+11, 3.79913127e+10}};

std::vector<double> res_exp =
    {4.633600473471703e+17};

double volume = (1.0e-27);
double dt = 1.0e-15;


template<typename RealType>
void test_expected_pair_number()
{
    for(int i = 0 ; i < E.size(); ++i){
        const auto exp = compute_expected_pair_number<RealType, unit_system::SI>(
            E[i][0], E[i][1], E[i][2], B[i][0], B[i][1], B[i][2],
            volume, dt);

        BOOST_CHECK_SMALL((exp-static_cast<RealType>(res_exp[i]))/static_cast<RealType>(res_exp[i]), tolerance<RealType>());
    }
}

BOOST_AUTO_TEST_CASE( picsar_schwinger_core_expected_pair_number )
{
    test_expected_pair_number <double>();
    test_expected_pair_number <float>();
}

/*

template<typename RealType>
constexpr void test_lambda_threshold_set_get()
{
    const auto t_lambda = static_cast<RealType>(lambda);

    auto rng_SI = wrapper_wrapper{765483};
    auto rng_norm_omega = wrapper_wrapper{11234};
    auto rng_norm_lambda = wrapper_wrapper{3409};

    auto engine_SI = schwinger_pair_engine<
        RealType, wrapper_wrapper , unit_system::SI>
        {std::move(rng_SI)};

    auto engine_norm_omega = schwinger_pair_engine<
        RealType, wrapper_wrapper , unit_system::norm_omega>
        {std::move(rng_norm_omega), t_lambda};

    auto engine_norm_lambda = schwinger_pair_engine<
        RealType, wrapper_wrapper , unit_system::norm_lambda>
        {std::move(rng_norm_lambda), t_lambda};

    BOOST_CHECK_EQUAL(engine_SI.get_lambda(), static_cast<RealType>(1.0));
    BOOST_CHECK_EQUAL(engine_norm_omega.get_lambda(), t_lambda);
    BOOST_CHECK_EQUAL(engine_norm_lambda.get_lambda(), t_lambda);

    const auto new_lambda = static_cast<RealType>(2.0*lambda);
    engine_SI.set_lambda(new_lambda);
    engine_norm_omega.set_lambda(new_lambda);
    engine_norm_lambda.set_lambda(new_lambda);

    BOOST_CHECK_EQUAL(engine_SI.get_lambda(), new_lambda);
    BOOST_CHECK_EQUAL(engine_norm_omega.get_lambda(), new_lambda);
    BOOST_CHECK_EQUAL(engine_norm_lambda.get_lambda(), new_lambda);

    BOOST_CHECK_EQUAL(engine_SI.get_poisson_gaussian_threshold(),
        schwinger_pair_engine<RealType, wrapper_wrapper , unit_system::SI>::poisson_gaussian_default_threshold);
    BOOST_CHECK_EQUAL(engine_norm_omega.get_poisson_gaussian_threshold(),
        schwinger_pair_engine<RealType, wrapper_wrapper , unit_system::norm_omega>
            ::poisson_gaussian_default_threshold);
    BOOST_CHECK_EQUAL(engine_norm_lambda.get_poisson_gaussian_threshold(),
        schwinger_pair_engine<RealType, wrapper_wrapper , unit_system::norm_lambda>
            ::poisson_gaussian_default_threshold);

    const auto new_threshold = static_cast<RealType>(50);
    engine_SI.set_poisson_gaussian_threshold(new_threshold);
    engine_norm_omega.set_poisson_gaussian_threshold(new_threshold);
    engine_norm_lambda.set_poisson_gaussian_threshold(new_threshold);

    BOOST_CHECK_EQUAL(engine_SI.get_poisson_gaussian_threshold(), new_threshold);
    BOOST_CHECK_EQUAL(engine_norm_omega.get_poisson_gaussian_threshold(), new_threshold);
    BOOST_CHECK_EQUAL(engine_norm_lambda.get_poisson_gaussian_threshold(), new_threshold);
}


BOOST_AUTO_TEST_CASE( picsar_schwinger_lambda_threshold_set_get )
{
    test_lambda_threshold_set_get <double>();
    test_lambda_threshold_set_get <float>();
}

*/
