//####### Test module for vec functions ####################################

//Define Module name
 #define BOOST_TEST_MODULE "math/poisson_distrib"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>

#include <array>
#include <utility>

#include "poisson_distrib.hpp"

using namespace picsar::multi_physics::math;

// ------------- Tests --------------


template<typename T>
void test_poisson()
{

    const auto lambda = static_cast<T>(5.0);


    const auto test_cases = std::array<std::pair<T, size_t>,16>{
        std::make_pair<T, size_t>(0.0, 0),
        std::make_pair<T, size_t>(.006, 0),
        std::make_pair<T, size_t>(.040, 1),
        std::make_pair<T, size_t>(.124, 2),
        std::make_pair<T, size_t>(.264, 3),
        std::make_pair<T, size_t>(.44, 4),
        std::make_pair<T, size_t>(.615,5),
        std::make_pair<T, size_t>(.76, 6),
        std::make_pair<T, size_t>(.86, 7),
        std::make_pair<T, size_t>(.93, 8),
        std::make_pair<T, size_t>(.96, 9),
        std::make_pair<T, size_t>(.98, 10),
        std::make_pair<T, size_t>(.99, 11),
        std::make_pair<T, size_t>(.997, 12),
        std::make_pair<T, size_t>(.999, 13),
        std::make_pair<T, size_t>(.9995, 14)};

    for (const auto& tc : test_cases){
        BOOST_CHECK_EQUAL(poisson_distrib(lambda, tc.first), tc.second);
    }
}

//Test empty constructor
BOOST_AUTO_TEST_CASE( picsar_poisson_distrib )
{
    test_poisson<double>();
    test_poisson<float>();
}
