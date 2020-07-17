//####### Test module for serialization ####################################

//Define Module name
 #define BOOST_TEST_MODULE "utils/serialization"

#include <vector>

//Include Boost unit tests library
#include <boost/test/unit_test.hpp>

#include "serialization.hpp"

using namespace picsar::multi_physics::utils::serialization;

// ------------- Tests --------------

// ***Test "put_in"

BOOST_AUTO_TEST_CASE( picsar_serialization_put_in)
{
    std::vector<char> raw_data;
    const int ii = -8291;
    const size_t ss = 9911;
    const double data_double = 1e12;
    const float data_float = -1e3;
    const char cc = '@';

    put_in(ii, raw_data);
    put_in(ss, raw_data);
    put_in(data_double, raw_data);
    put_in(data_float, raw_data);
    put_in(cc, raw_data);

    char* ptr_raw_data = raw_data.data();
    size_t index = 0;

    BOOST_CHECK_EQUAL(
        *reinterpret_cast<int*>(ptr_raw_data+index),
        ii);
    index += sizeof(int);

    BOOST_CHECK_EQUAL(
        *reinterpret_cast<size_t*>(ptr_raw_data+index),
        ss);
    index += sizeof(size_t);

    BOOST_CHECK_EQUAL(
        *reinterpret_cast<double*>(ptr_raw_data+index),
        data_double);
    index += sizeof(double);

    BOOST_CHECK_EQUAL(
        *reinterpret_cast<float*>(ptr_raw_data+index),
        data_float);
    index += sizeof(float);

    BOOST_CHECK_EQUAL(
        *reinterpret_cast<char*>(ptr_raw_data+index),
        cc);
}

// *******************************

// ***Test "get_out"

BOOST_AUTO_TEST_CASE( picsar_serialization_get_out)
{
    std::vector<char> raw_data;

    const int ii = -8291;
    const size_t ss = 9911;
    const double data_double = 1e12;
    const float data_float = -1e3;
    const char cc = '@';

    put_in(ii, raw_data);
    put_in(ss, raw_data);
    put_in(data_double, raw_data);
    put_in(data_float, raw_data);
    put_in(cc, raw_data);

    auto it = raw_data.begin();

    BOOST_CHECK_EQUAL(get_out<int>(it), ii);
    BOOST_CHECK_EQUAL(get_out<size_t>(it), ss);
    BOOST_CHECK_EQUAL(get_out<double>(it), data_double);
    BOOST_CHECK_EQUAL(get_out<float>(it), data_float);
    BOOST_CHECK_EQUAL(get_out<char>(it), cc);
}

// *******************************

// ***Test "get_out" (version for multiple variables)

BOOST_AUTO_TEST_CASE( picsar_serialization_get_n_out)
{
    std::vector<char> raw_data;

    auto doublevec = std::vector<double>{
        1.0, 2,0, 3.0, 4.0, 5.0, 6.0};

    for (auto dd : doublevec)
        put_in(dd, raw_data);

    auto it = raw_data.begin();
    auto dvec1 = get_n_out<double>(it, 3);
    auto dvec2 = get_n_out<double>(it, 3);

    for(int i = 0; i < 3; ++i){
        BOOST_CHECK_EQUAL(dvec1[i], doublevec[i]);
        BOOST_CHECK_EQUAL(dvec2[i], doublevec[i+3]);
    }
}

// *******************************
