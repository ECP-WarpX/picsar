//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "containers/picsar_tables"

//Will automatically define a main for this test
 #define BOOST_TEST_DYN_LINK

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>
#include <algorithm>
#include <array>
#include <iostream>

#include "picsar_tables.hpp"

using namespace picsar::multi_physics::containers;

//Tolerance for double precision calculations
const double double_tolerance = 1.0e-12;

//Tolerance for single precision calculations
const float float_tolerance = 1.0e-5;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

double linear_function(double x)
{
    const double m = -3.0;
    const double q = 1.7;
    return m*x + q;
}

const double xmin = -10.0;
const double xmax = 10.0;
const int size = 100;

equispaced_1d_table<double, std::vector<double> > make_1d_table()
{
    std::vector<double> data(size);
    std::generate(data.begin(), data.end(),
        [&, n = 0]() mutable {
            double x = xmin+((xmax-xmin)*(n++))/(size-1);
            return linear_function(x);
        });

    return equispaced_1d_table<double, std::vector<double>>(xmin, xmax,data);
}

// ------------- Tests --------------

void check_table_1d(
    const equispaced_1d_table<double, std::vector<double> >& tab)
{
    BOOST_CHECK_EQUAL(0, 0);

    const auto rxmin = tab.get_x_min();
    BOOST_CHECK_EQUAL(rxmin, xmin);
    const auto rxmax = tab.get_x_max();
    BOOST_CHECK_EQUAL(rxmax, xmax);
    const auto rhowmany =  tab.get_how_many_x();
    BOOST_CHECK_EQUAL(rhowmany, size);

    const auto x0 = tab.x_coord(0);
    const auto x1 = tab.x_coord(size/2);
    const auto x2 = tab.x_coord(size-1);
    const auto x1exp = (size/2)*(xmax - xmin)/(size-1) + xmin;
    BOOST_CHECK_SMALL(fabs((x0-xmin)/xmin), tolerance<double>());
    BOOST_CHECK_SMALL(fabs(x1 - x1exp)/x1exp, tolerance<double>());
    BOOST_CHECK_SMALL(fabs((x2 - xmax)/xmax), tolerance<double>());

    const auto x3 = 0.732*(x1-x0) + x0;
    const auto x4 = 0.118*(x2-x1) + x1;
    const auto xarr = std::array<double, 5>{x0,x1,x2,x3,x4};
    for (const auto& xx : xarr)
    {
        const auto val = tab.interp(xx);
        const auto expected = linear_function(xx);
        BOOST_CHECK_SMALL(fabs((val - expected)/expected), tolerance<double>());
    }
}

// ***Test equispaced_1d_table constructor and getters
BOOST_AUTO_TEST_CASE( picsar_equispaced_1d_table_constructor_getters)
{
    auto tab_1d = make_1d_table();
    const auto const_tab_1d = make_1d_table();
    auto copy_tab_1d = tab_1d;

    check_table_1d(tab_1d);
    check_table_1d(const_tab_1d);
    check_table_1d(copy_tab_1d);
}

// ***Test equispaced_1d_table setter
BOOST_AUTO_TEST_CASE( picsar_equispaced_1d_table_constructor_setter)
{
    auto tab_1d = make_1d_table();
    const auto val = 1000.0;
    const int where = 10;
    tab_1d.set_val(where, val);
    const auto x10 = tab_1d.x_coord(where);
    const auto res =tab_1d.interp(x10);

    BOOST_CHECK_SMALL(fabs((res - val)/val), tolerance<double>());

}

// *******************************
