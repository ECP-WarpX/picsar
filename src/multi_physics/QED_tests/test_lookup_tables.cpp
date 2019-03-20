//####### Test module for lookup tables ######################################

//Define Module name
 #define BOOST_TEST_MODULE "lookup tables"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

#include <cmath>
#include <array>
#include <vector>
#include <functional>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//Units choice. Not relevant here, but avoids compile-time warning
#define PXRMP_USE_SI_UNITS

#include "lookup_tables.hpp"

using namespace picsar::multi_physics;

// ------------- Tests --------------

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-10;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-4;

//Templated tolerance
template <typename T>
T tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//________________________1D tables__________________________________

//Test lookup_1d constructor generic
template<typename T>
void lookup1d_constructor()
{
    std::vector<T> coords{0,1,2,3,4};
    std::vector<T> vals;

    std::transform(coords.begin(), coords.end(), std::back_inserter(vals),
        [](T val){return static_cast<T>(2.0*val);});

    lookup_1d<T> l1d{coords, vals,
        lookup_1d<T>::linear_equispaced_interpolation};

    lookup_1d<T> l1dbis{l1d};

    lookup_1d<T> l1dnull;

    lookup_1d<T> l1copy = l1d;



    for(size_t i = 0; i < coords.size(); ++i){
        BOOST_CHECK_EQUAL(l1d.get_coords()[i], static_cast<T>(coords[i]));
        BOOST_CHECK_EQUAL(l1dbis.get_coords()[i], static_cast<T>(coords[i]));
        BOOST_CHECK_EQUAL(l1copy.get_coords()[i], static_cast<T>(coords[i]));
    }

    BOOST_CHECK_EQUAL(l1d.is_init(),true);
    BOOST_CHECK_EQUAL(l1dbis.is_init(),true);
    BOOST_CHECK_EQUAL(l1dnull.is_init(),false);
    BOOST_CHECK_EQUAL(l1copy.is_init(),true);
}

//Test lookup_1d constructor in double precision
BOOST_AUTO_TEST_CASE( lookup1d_constructor_double_1 )
{
    lookup1d_constructor<double>();
}

//Test lookup_1d constructor in single precision
BOOST_AUTO_TEST_CASE( lookup1d_constructor_single_1 )
{
    lookup1d_constructor<float>();
}

//Test lookup_1d linear equispaced interpolator generic
template<typename T>
void lookup1d_linear_equi_interp()
{
    std::vector<T> coords{0,1,2,3,4};
    std::vector<T> vals;

    std::transform(coords.begin(), coords.end(), std::back_inserter(vals),
        [](T val){return static_cast<T>(2.0*val);});

    lookup_1d<T> l1d{coords, vals,
        lookup_1d<T>::linear_equispaced_interpolation};

    std::vector<T> where{0.0001, 0.5, 1.0, 1.5, 2.0, 3.5, 3.999};

    for(auto ww: where){
        T res = l1d.interp(ww);
        T exp = static_cast<T>(2.0*ww);
        BOOST_CHECK_SMALL((res-exp)/(exp), tolerance<T>());
    }

    BOOST_CHECK_SMALL(l1d.interp(static_cast<T>(0.0)), tolerance<T>());

}

//Test lookup_1d constructor in double precision
BOOST_AUTO_TEST_CASE( lookup1d_linear_interp_equi_double_1 )
{
    lookup1d_linear_equi_interp<double>();
}

//Test lookup_1d constructor in single precision
BOOST_AUTO_TEST_CASE( lookup1d_linear_interp_equi_single_1 )
{
    lookup1d_linear_equi_interp<float>();
}


//Test lookup_1d linear interpolator generic
template<typename T>
void lookup1d_linear_interp()
{
    std::vector<T> coords{0,1,2,3,4};
    std::vector<T> vals;

    std::transform(coords.begin(), coords.end(), std::back_inserter(vals),
        [](T val){return static_cast<T>(2.0*val);});

    lookup_1d<T> l1d{coords, vals,
        lookup_1d<T>::linear_interpolation};

    std::vector<T> where{0.0001, 0.5, 1.0, 1.5, 2.0, 3.5, 3.999};

    for(auto ww: where){
        T res = l1d.interp(ww);
        T exp = static_cast<T>(2.0*ww);
        BOOST_CHECK_SMALL((res-exp)/(exp), tolerance<T>());
    }

    BOOST_CHECK_SMALL(l1d.interp(static_cast<T>(0.0)), tolerance<T>());

}

//Test lookup_1d constructor in double precision
BOOST_AUTO_TEST_CASE( lookup1d_linear_interp_double_1 )
{
    lookup1d_linear_interp<double>();
}

//Test lookup_1d constructor in single precision
BOOST_AUTO_TEST_CASE( lookup1d_linear_interp_single_1 )
{
    lookup1d_linear_interp<float>();
}


//________________________2D tables__________________________________

//Test lookup_2d constructor generic
template<typename T>
void lookup2d_constructor()
{
    std::vector<T> c1{-3,-2,-1,0,1,2,3};
    std::vector<T> c2{-3,-2,-1,0,1,2,3};
    std::vector<T> vals{};

    for (auto cc1: c1){
        for(auto cc2: c2){
            vals.push_back(cc1+cc2);
        }
    }

    lookup_2d<T> l2d{std::array<std::vector<T>,2>{c1,c2}, vals,
        lookup_2d<T>::linear_interpolation, lookup_2d<T>::row_major};

    lookup_2d<T> l2dbis{l2d};

    lookup_2d<T> l2dnull;

    lookup_2d<T> l2copy = l2d;

    for(size_t i = 0; i < c1.size(); ++i){
        BOOST_CHECK_EQUAL(l2d.get_coords()[0][i], static_cast<T>(c1[i]));
        BOOST_CHECK_EQUAL(l2dbis.get_coords()[0][i], static_cast<T>(c1[i]));
        BOOST_CHECK_EQUAL(l2copy.get_coords()[0][i], static_cast<T>(c1[i]));
    }

    for(size_t i = 0; i < c2.size(); ++i){
        BOOST_CHECK_EQUAL(l2d.get_coords()[1][i], static_cast<T>(c2[i]));
        BOOST_CHECK_EQUAL(l2dbis.get_coords()[1][i], static_cast<T>(c2[i]));
        BOOST_CHECK_EQUAL(l2copy.get_coords()[1][i], static_cast<T>(c2[i]));
    }

    BOOST_CHECK_EQUAL(l2d.is_init(),true);
    BOOST_CHECK_EQUAL(l2dbis.is_init(),true);
    BOOST_CHECK_EQUAL(l2dnull.is_init(),false);
    BOOST_CHECK_EQUAL(l2copy.is_init(),true);
}

//Test lookup_1d constructor in double precision
BOOST_AUTO_TEST_CASE( lookup2d_constructor_double_1 )
{
    lookup2d_constructor<double>();
}

//Test lookup_1d constructor in single precision
BOOST_AUTO_TEST_CASE( lookup2d_constructor_single_1 )
{
    lookup2d_constructor<float>();
}

//Test lookup_2d linear interpolator generic
template<typename T>
void lookup2d_linear_interpolator()
{
    std::vector<T> c1{-3,-2,-1,0,1,2,3};
    std::vector<T> c2{-3,-2,-1,0,1,2,3};
    std::vector<T> vals{};

    for (auto cc1: c1){
        for(auto cc2: c2){
            vals.push_back(cc1+cc2);
        }
    }

    lookup_2d<T> l2d{std::array<std::vector<T>,2>{c1,c2}, vals,
        lookup_2d<T>::linear_interpolation, lookup_2d<T>::row_major};

    std::vector<T> c1_test{-2.9, -2.1, 1.2, 0.0, 2.99};
    std::vector<T> c2_test{-1.9, -2.7, 2.1, 1.1, 0.7};

    for(auto cc1: c1_test){
        for(auto cc2: c1_test){
            T res = l2d.interp(cc1, cc2);
            T exp = static_cast<T>(cc1+cc2);
            BOOST_CHECK_SMALL((res-exp)/(exp), tolerance<T>());
        }
    }
    BOOST_CHECK_EQUAL(1,1);
}

//Test lookup_2d linear_interpolator in double precision
BOOST_AUTO_TEST_CASE( lookup2d_linear_interpolator_double_1 )
{
    lookup2d_linear_interpolator<double>();
}

//Test lookup_2d linear_interpolator in single precision
BOOST_AUTO_TEST_CASE( lookup2d_linear_interpolator_single_1 )
{
    lookup2d_linear_interpolator<float>();
}
