//####### Test module for picsar_tables ####################################

//Define Module name
 #define BOOST_TEST_MODULE "containers/picsar_tables"

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <picsar_qed/containers/picsar_tables.hpp>
#include <picsar_qed/containers/picsar_span.hpp>
#include <picsar_qed/utils/serialization.hpp>

#include <vector>
#include <algorithm>
#include <array>

using namespace picsar::multi_physics::containers;

namespace pxr_ser = picsar::multi_physics::utils::serialization;

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

double linear_function(double x, double y)
{
    const double a = 2.2;
    const double b = -1.7;
    const double c = -3.1;
    return a*x + b*y + c;
}

const double xmin = -10.0;
const double xmax = 10.0;
const double ymin = -8.0;
const double ymax = 8.0;
const int xsize = 100;
const int ysize = 100;

const int gxsize = 105;
const int gysize = 105;
const double gxmin = 1e-12;
const double gxmax = 11.0;
const double gymin = 1e-10;
const double gymax = 9.0;
const double gxswitch = 8.0;
const double gyswitch = 6.0;
const int gxfirst = 80;
const int gyfirst = 80;

class Functor
{
    public:

    Functor(){}

    Functor(const int zsize, const int zfirst,
        const double zmin, const double zmax, const double zswitch):
        m_zsize{zsize}, m_zfirst{zfirst},
        m_zmin{zmin}, m_zmax{zmax}, m_zswitch{zswitch}
    {
        m_logzmin = std::log(m_zmin);
        m_logzswitch = std::log(m_zswitch);
    }

    double operator() (const int i) const
    {
        if (i < m_zfirst){
            return std::exp(m_logzmin + i*(m_logzswitch-m_logzmin)/(m_zfirst-1));
        }
        else{
            return m_zswitch+((i+1)-m_zfirst)*(m_zmax-m_zswitch)/(m_zsize-m_zfirst);
        }
    }

    bool operator== (const Functor& rhs) const
    {
        return
            (this->m_zsize == rhs.m_zsize) &&
            (this->m_zfirst == rhs.m_zfirst) &&
            (this->m_zmin == rhs.m_zmin) &&
            (this->m_logzmin == rhs.m_logzmin) &&
            (this->m_zmax == rhs.m_zmax) &&
            (this->m_zswitch == rhs.m_zswitch) &&
            (this->m_logzswitch == rhs.m_logzswitch);
    }

    std::vector<char> serialize() const
    {
        std::vector<char> raw;

        pxr_ser::put_in(m_zsize, raw);
        pxr_ser::put_in(m_zfirst, raw);
        pxr_ser::put_in(m_zmin, raw);
        pxr_ser::put_in(m_logzmin, raw);
        pxr_ser::put_in(m_zmax, raw);
        pxr_ser::put_in(m_zswitch, raw);
        pxr_ser::put_in(m_logzswitch, raw);

        return raw;
    }

    template<typename Iter>
    static Functor deserialize(Iter& it)
    {
        auto functor = Functor{};

        functor.m_zsize = pxr_ser::get_out<int>(it);
        functor.m_zfirst = pxr_ser::get_out<int>(it);
        functor.m_zmin = pxr_ser::get_out<double>(it);
        functor.m_logzmin = pxr_ser::get_out<double>(it);
        functor.m_zmax = pxr_ser::get_out<double>(it);
        functor.m_zswitch = pxr_ser::get_out<double>(it);
        functor.m_logzswitch = pxr_ser::get_out<double>(it);

        return functor;
    }

    private:

    int m_zsize;
    int m_zfirst;
    double m_zmin;
    double m_logzmin;
    double m_zmax;
    double m_zswitch;
    double m_logzswitch;
};

class IFunctor
{
    public:

    IFunctor(){}

    IFunctor(const int zsize, const int zfirst,
        const double zmin, const double zmax, const double zswitch):
        m_zsize{zsize}, m_zfirst{zfirst},
        m_zmin{zmin}, m_zmax{zmax}, m_zswitch{zswitch}
    {
        m_logzmin = std::log(m_zmin);
        m_logzswitch = std::log(m_zswitch);
    }

    int operator() (const double z) const
    {
        if (z < m_zswitch){
            const auto logz = std::log(z);
            return static_cast<int>(
                std::round( (m_zfirst -1)*(logz-m_logzmin)/(m_logzswitch - m_logzmin)));
        }
        else{
            return static_cast<int>(
                std::floor(m_zfirst + (m_zsize-m_zfirst - 1)*(z-m_zswitch)/(m_zmax - m_zswitch)));
        }
    }

    bool operator== (const IFunctor& rhs) const
    {
        return
            (this->m_zsize == rhs.m_zsize) &&
            (this->m_zfirst == rhs.m_zfirst) &&
            (this->m_zmin == rhs.m_zmin) &&
            (this->m_logzmin == rhs.m_logzmin) &&
            (this->m_zmax == rhs.m_zmax) &&
            (this->m_zswitch == rhs.m_zswitch) &&
            (this->m_logzswitch == rhs.m_logzswitch);
    }

    std::vector<char> serialize() const
    {
        std::vector<char> raw;

        pxr_ser::put_in(m_zsize, raw);
        pxr_ser::put_in(m_zfirst, raw);
        pxr_ser::put_in(m_zmin, raw);
        pxr_ser::put_in(m_logzmin, raw);
        pxr_ser::put_in(m_zmax, raw);
        pxr_ser::put_in(m_zswitch, raw);
        pxr_ser::put_in(m_logzswitch, raw);

        return raw;
    }

    template<typename Iter>
    static IFunctor deserialize(Iter& it)
    {
        auto functor = IFunctor{};

        functor.m_zsize = pxr_ser::get_out<int>(it);
        functor.m_zfirst = pxr_ser::get_out<int>(it);
        functor.m_zmin = pxr_ser::get_out<double>(it);
        functor.m_logzmin = pxr_ser::get_out<double>(it);
        functor.m_zmax = pxr_ser::get_out<double>(it);
        functor.m_zswitch = pxr_ser::get_out<double>(it);
        functor.m_logzswitch = pxr_ser::get_out<double>(it);

        return functor;
    }

    private:

    int m_zsize;
    int m_zfirst;
    double m_zmin;
    double m_logzmin;
    double m_zmax;
    double m_zswitch;
    double m_logzswitch;
};

const auto x_functor = Functor{gxsize, gxfirst, gxmin, gxmax, gxswitch};
const auto y_functor = Functor{gysize, gyfirst, gymin, gymax, gyswitch};
const auto ix_functor = IFunctor{gxsize, gxfirst, gxmin, gxmax, gxswitch};
const auto iy_functor = IFunctor{gysize, gyfirst, gymin, gymax, gyswitch};

equispaced_1d_table<double, std::vector<double> > make_1d_table()
{
    std::vector<double> data(xsize);
    std::generate(data.begin(), data.end(),
        [&, n = 0]() mutable {
            double x = xmin+((xmax-xmin)*(n++))/(xsize-1);
            return linear_function(x);
        });

    return equispaced_1d_table<double, std::vector<double>>(xmin, xmax,data);
}

equispaced_2d_table<double, std::vector<double> > make_2d_table()
{
    std::vector<double> data(xsize*ysize);
    for (int i = 0; i < xsize; ++i){
        for (int j = 0; j < ysize; ++j){
            double x = xmin+i*(xmax-xmin)/(xsize-1);
            double y = ymin+j*(ymax-ymin)/(ysize-1);
            data[i*ysize+j] = linear_function(x, y);
        }
    }

    return equispaced_2d_table<double, std::vector<double>>(
        xmin, xmax, ymin, ymax, xsize, ysize, data);
}

using Generic2DTableType = generic_2d_table<
        double, std::vector<double>,
        Functor, Functor,
        IFunctor, IFunctor>;

auto make_generic_2d_table()
{
    auto tab = Generic2DTableType(
            gxsize, gysize,
            std::vector<double>(gxsize*gysize),
            x_functor, y_functor,
            ix_functor, iy_functor);

    const auto coords = tab.get_all_coordinates();
    auto vals = std::vector<double>(coords.size());

    int i = 0;
    for (const auto cc : coords){
        tab.set_val(i++, linear_function(cc[0], cc[1]));
    }

    return tab;
}


// ------------- Tests --------------

void check_table_1d(
    const equispaced_1d_table<double, std::vector<double> >& tab)
{
    const auto rxmin = tab.get_x_min();
    BOOST_CHECK_EQUAL(rxmin, xmin);
    const auto rxmax = tab.get_x_max();
    BOOST_CHECK_EQUAL(rxmax, xmax);
    const auto rxsize = tab.get_x_size();
    BOOST_CHECK_EQUAL(rxsize, xmax-xmin);
    const auto rhowmany =  tab.get_how_many_x();
    BOOST_CHECK_EQUAL(rhowmany, xsize);
    const auto rdx = tab.get_dx();
    BOOST_CHECK_EQUAL(rdx, (xmax-xmin)/(xsize-1));

    const auto first_val = tab.get_values_reference()[0];
    BOOST_CHECK_EQUAL(first_val, linear_function(xmin));

    const auto x0 = tab.get_x_coord(0);
    const auto x1 = tab.get_x_coord(xsize/2);
    const auto x2 = tab.get_x_coord(xsize-1);
    const auto x1exp = (xsize/2)*(xmax - xmin)/(xsize-1) + xmin;
    BOOST_CHECK_SMALL((x0-xmin)/xmin, tolerance<double>());
    BOOST_CHECK_SMALL((x1 - x1exp)/x1exp, tolerance<double>());
    BOOST_CHECK_SMALL((x2 - xmax)/xmax, tolerance<double>());

    const auto x3 = 0.732*(x1-x0) + x0;
    const auto x4 = 0.118*(x2-x1) + x1;

    const auto v0 = tab.get_val(0);
    const auto v1 = tab.get_val(xsize/2);
    const auto v2 = tab.get_val(xsize-1);
    BOOST_CHECK_EQUAL(v0, linear_function(x0));
    BOOST_CHECK_EQUAL(v1, linear_function(x1));
    BOOST_CHECK_EQUAL(v2, linear_function(x2));

    const auto xarr = std::array<double, 5>{x0,x1,x2,x3,x4};
    for (const auto& xx : xarr)
    {
        const auto val = tab.interp(xx);
        const auto expected = linear_function(xx);
        BOOST_CHECK_SMALL((val - expected)/expected, tolerance<double>());
    }

    const auto all_coords = tab.get_all_coordinates();
    BOOST_CHECK_EQUAL(all_coords.size(), xsize);
    for (int i = 0; i < xsize; ++i){
        BOOST_CHECK_EQUAL(all_coords[i], tab.get_x_coord(i));
    }
}

void check_table_2d(
    const equispaced_2d_table<double, std::vector<double> >& tab)
{
    const auto rxmin = tab.get_x_min();
    BOOST_CHECK_EQUAL(rxmin, xmin);
    const auto rxmax = tab.get_x_max();
    BOOST_CHECK_EQUAL(rxmax, xmax);
    const auto rhowmany_x =  tab.get_how_many_x();
    BOOST_CHECK_EQUAL(rhowmany_x, xsize);
    const auto rxsize =  tab.get_x_size();
    BOOST_CHECK_EQUAL(rxsize, xmax-xmin);
    const auto rdx =  tab.get_dx();
    BOOST_CHECK_EQUAL(rdx, (xmax-xmin)/(rhowmany_x-1));
    const auto rymin = tab.get_y_min();
    BOOST_CHECK_EQUAL(rymin, ymin);
    const auto rymax = tab.get_y_max();
    BOOST_CHECK_EQUAL(rymax, ymax);
    const auto rhowmany_y =  tab.get_how_many_y();
    BOOST_CHECK_EQUAL(rhowmany_y, ysize);
    const auto rysize =  tab.get_y_size();
    BOOST_CHECK_EQUAL(rysize, ymax-ymin);
    const auto rdy =  tab.get_dy();
    BOOST_CHECK_EQUAL(rdy, (ymax-ymin)/(rhowmany_y-1));

    const auto first_val = tab.get_values_reference()[0];
    BOOST_CHECK_EQUAL(first_val, linear_function(xmin, ymin));

    const auto x0 = tab.get_x_coord(0);
    const auto x1 = tab.get_x_coord(xsize/2);
    const auto x2 = tab.get_x_coord(xsize-1);
    const auto x1exp = (xsize/2)*(xmax - xmin)/(xsize-1) + xmin;
    BOOST_CHECK_SMALL((x0-xmin)/xmin, tolerance<double>());
    BOOST_CHECK_SMALL((x1 - x1exp)/x1exp, tolerance<double>());
    BOOST_CHECK_SMALL((x2 - xmax)/xmax, tolerance<double>());

    const auto x3 = 0.732*(x1-x0) + x0;
    const auto x4 = 0.118*(x2-x1) + x1;
    const auto xarr = std::array<double, 5>{x0,x1,x2,x3,x4};

    const auto y0 = tab.get_y_coord(0);
    const auto y1 = tab.get_y_coord(ysize/2);
    const auto y2 = tab.get_y_coord(ysize-1);
    const auto y1exp = (ysize/2)*(ymax - ymin)/(ysize-1) + ymin;
    BOOST_CHECK_SMALL((y0-ymin)/ymin, tolerance<double>());
    BOOST_CHECK_SMALL((y1 - y1exp)/y1exp, tolerance<double>());
    BOOST_CHECK_SMALL((y2 - ymax)/ymax, tolerance<double>());

    const auto v00 = tab.get_val(0, 0);
    const auto v10 = tab.get_val(xsize/2, 0);
    const auto v20 = tab.get_val(xsize-1, 0);
    const auto v01 = tab.get_val(0, ysize/2);
    const auto v11 = tab.get_val(xsize/2, ysize/2);
    const auto v21 = tab.get_val(xsize-1, ysize/2);
    const auto v02 = tab.get_val(0, ysize-1);
    const auto v12 = tab.get_val(xsize/2, ysize-1);
    const auto v22 = tab.get_val(xsize-1, ysize-1);
    BOOST_CHECK_EQUAL(v00, linear_function(x0,y0));
    BOOST_CHECK_EQUAL(v10, linear_function(x1,y0));
    BOOST_CHECK_EQUAL(v20, linear_function(x2,y0));
    BOOST_CHECK_EQUAL(v01, linear_function(x0,y1));
    BOOST_CHECK_EQUAL(v11, linear_function(x1,y1));
    BOOST_CHECK_EQUAL(v21, linear_function(x2,y1));
    BOOST_CHECK_EQUAL(v02, linear_function(x0,y2));
    BOOST_CHECK_EQUAL(v12, linear_function(x1,y2));
    BOOST_CHECK_EQUAL(v22, linear_function(x2,y2));

    const auto y3 = 0.8569*(y1-y0) + y0;
    const auto y4 = 0.3467*(y2-y1) + y1;
    const auto yarr = std::array<double, 5>{y0,y1,y2,y3,y4};

    for (const auto& xx : xarr)
    {
        for (const auto& yy : yarr)
        {
            const auto val = tab.interp(xx, yy);
            const auto expected = linear_function(xx, yy);
            BOOST_CHECK_SMALL((val - expected)/expected, tolerance<double>());
        }
    }

    const auto all_coords = tab.get_all_coordinates();
    int count = 0;
    BOOST_CHECK_EQUAL(all_coords.size(), xsize*ysize);
    for (int i = 0; i < xsize; ++i){
        for (int j = 0; j < ysize; ++j){
            auto cc = all_coords[count++];
            BOOST_CHECK_EQUAL(cc[0], tab.get_x_coord(i));
            BOOST_CHECK_EQUAL(cc[1], tab.get_y_coord(j));
        }
    }
}

void check_table_2d_interp_one_coord(
        const equispaced_2d_table<double, std::vector<double> >& tab)
{
    const auto x0 = tab.get_x_coord(0);
    const auto x1 = tab.get_x_coord(xsize/2);
    const auto x2 = tab.get_x_coord(xsize-1);
    const auto x3 = 0.732*(x1-x0) + x0;
    const auto x4 = 0.118*(x2-x1) + x1;
    const auto xarr = std::array<double, 5>{x0,x1,x2,x3,x4};
    const auto jarr = std::array<int, 5>{0,1,5,27,ysize-1};

    for (auto xx : xarr){
        for (auto jj : jarr)
        {
            const auto yy = tab.get_y_coord(jj);
            const auto res = tab.interp_first_coord(xx, jj);
            const auto exp = linear_function(xx, yy);
            BOOST_CHECK_SMALL((res - exp)/exp, tolerance<double>());
        }
    }


    const auto y0 = tab.get_y_coord(0);
    const auto y1 = tab.get_y_coord(ysize/2);
    const auto y2 = tab.get_y_coord(ysize-1);
    const auto y3 = 0.8569*(y1-y0) + y0;
    const auto y4 = 0.3467*(y2-y1) + y1;
    const auto yarr = std::array<double, 5>{y0,y1,y2,y3,y4};
    const auto iarr = std::array<int, 5>{0,1,5,27,xsize-1};

    for (auto yy : yarr){
        for (auto ii : iarr)
        {
            const auto xx = tab.get_x_coord(ii);
            const auto res = tab.interp_second_coord(ii, yy);
            const auto exp = linear_function(xx, yy);
            BOOST_CHECK_SMALL((res - exp)/exp, tolerance<double>());
        }
    }
}

void check_generic_table_2d(
    const Generic2DTableType& tab)
{
    const auto rxmin = tab.get_x_min();
    BOOST_CHECK_SMALL((rxmin - gxmin)/gxmin, tolerance<double>());
    const auto rxmax = tab.get_x_max();
    BOOST_CHECK_SMALL((rxmax - gxmax)/gxmax, tolerance<double>());
    const auto rhowmany_x =  tab.get_how_many_x();
    BOOST_CHECK_EQUAL(rhowmany_x, gxsize);
    const auto rxsize =  tab.get_x_size();
    BOOST_CHECK_EQUAL(rxsize, gxmax-gxmin);
    const auto rymin = tab.get_y_min();
    BOOST_CHECK_SMALL((rymin - gymin)/gymin, tolerance<double>());
    const auto rymax = tab.get_y_max();
    BOOST_CHECK_SMALL((rymax - gymax)/gymax, tolerance<double>());
    const auto rhowmany_y =  tab.get_how_many_y();
    BOOST_CHECK_EQUAL(rhowmany_y, gysize);
    const auto rysize =  tab.get_y_size();
    BOOST_CHECK_EQUAL(rysize, gymax-gymin);
}

void check_generic_table_2d_interp_one_coord(
        const Generic2DTableType& tab)
{
    const auto x0 = tab.get_x_coord(0);
    const auto x1 = tab.get_x_coord(xsize/4);
    const auto x2 = tab.get_x_coord(xsize-1);
    const auto x3 = 0.732*(x1-x0) + x0;
    const auto x4 = 0.118*(x2-x1) + x1;
    const auto xarr = std::array<double, 5>{x0,x1,x2,x3,x4};
    const auto jarr = std::array<int, 5>{0,1,5,27,ysize-1};

    for (auto xx : xarr){
        for (auto jj : jarr)
        {
            const auto yy = tab.get_y_coord(jj);
            const auto res = tab.interp_first_coord(xx, jj);
            const auto exp = linear_function(xx, yy);
            BOOST_CHECK_SMALL((res - exp)/exp, tolerance<double>());
        }
    }


    const auto y0 = tab.get_y_coord(0);
    const auto y1 = tab.get_y_coord(ysize/4);
    const auto y2 = tab.get_y_coord(ysize-1);
    const auto y3 = 0.8569*(y1-y0) + y0;
    const auto y4 = 0.3467*(y2-y1) + y1;
    const auto yarr = std::array<double, 5>{y0,y1,y2,y3,y4};
    const auto iarr = std::array<int, 5>{0,1,5,27,xsize-1};

    for (auto yy : yarr){
        for (auto ii : iarr)
        {
            const auto xx = tab.get_x_coord(ii);
            const auto res = tab.interp_second_coord(ii, yy);
            const auto exp = linear_function(xx, yy);
            BOOST_CHECK_SMALL((res - exp)/exp, tolerance<double>());
        }
    }
}

void check_generic_table_2d_interp(int i, int j, double val)
{
    auto gtab_2d = make_generic_2d_table();
    gtab_2d.set_val(i,j,val);
    const auto xx = gtab_2d.get_x_coord(i);
    const auto yy = gtab_2d.get_y_coord(j);
    const auto res = gtab_2d.interp(xx,yy);

    BOOST_CHECK_SMALL((res - val)/val, tolerance<double>());
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
BOOST_AUTO_TEST_CASE( picsar_equispaced_1d_table_setter)
{
    auto tab_1d = make_1d_table();
    const auto val = 1000.0;
    const int where = 10;
    tab_1d.set_val(where, val);
    const auto x10 = tab_1d.get_x_coord(where);
    const auto res =tab_1d.interp(x10);

    BOOST_CHECK_SMALL((res - val)/val, tolerance<double>());
}

// ***Test equispaced_1d_table set all
BOOST_AUTO_TEST_CASE( picsar_equispaced_1d_table_setall)
{
    auto tab_1d = make_1d_table();
    const auto val = 42.0;
    const auto vec = std::vector<double>(tab_1d.get_how_many_x(), val);

    tab_1d.set_all_vals(vec);

    const auto x0 = tab_1d.get_x_min();
    const auto x1 = tab_1d.get_x_max();
    const auto x2 = 0.376*(x1-x0) + x0;
    const auto x3 = 0.688*(x1-x0) + x0;

    const auto r0 =tab_1d.interp(x0);
    const auto r1 =tab_1d.interp(x1);
    const auto r2 =tab_1d.interp(x2);
    const auto r3 =tab_1d.interp(x3);

    BOOST_CHECK_SMALL((r0 - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((r1 - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((r2 - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((r3 - val)/val, tolerance<double>());
}

// ***Test equispaced_1d_table equality
BOOST_AUTO_TEST_CASE( picsar_equispaced_1d_table_equality)
{
    auto tab_1d = make_1d_table();
    auto tab_1d_2 = make_1d_table();
    BOOST_CHECK_EQUAL(tab_1d == tab_1d_2, true);

    tab_1d.set_val(1, 3.14);
    BOOST_CHECK_EQUAL(tab_1d == tab_1d_2, false);
}

// ***Test equispaced_1d_table serialization
BOOST_AUTO_TEST_CASE( picsar_equispaced_1d_table_serialization)
{
    auto tab_1d = make_1d_table();
    auto raw_data = tab_1d.serialize();
    auto tab_1d_2 = equispaced_1d_table<double, std::vector<double>>{raw_data};

    BOOST_CHECK_EQUAL(tab_1d.get_x_min(), tab_1d_2.get_x_min());
    BOOST_CHECK_EQUAL(tab_1d.get_x_max(), tab_1d_2.get_x_max());
    BOOST_CHECK_EQUAL(tab_1d.get_x_size(), tab_1d_2.get_x_size());
    BOOST_CHECK_EQUAL(tab_1d.get_how_many_x(), tab_1d_2.get_how_many_x());

    for(int i = 0; i < tab_1d.get_how_many_x(); ++i){
        BOOST_CHECK_EQUAL(tab_1d.get_val(i), tab_1d_2.get_val(i));
    }
}

// ***Test equispaced_2d_table constructor and getters
BOOST_AUTO_TEST_CASE( picsar_equispaced_2d_table_constructor_getters)
{
    auto tab_2d = make_2d_table();
    const auto const_tab_2d = make_2d_table();
    auto copy_tab_2d = tab_2d;

    check_table_2d(tab_2d);
    check_table_2d(const_tab_2d);
    check_table_2d(copy_tab_2d);
}

// ***Test equispaced_2d_table single coordinate interpolation
BOOST_AUTO_TEST_CASE( picsar_equispaced_2d_table_interp_one_coord)
{
    const auto const_tab_2d = make_2d_table();

    check_table_2d_interp_one_coord(const_tab_2d);
}

// ***Test equispaced_2d_table setter
BOOST_AUTO_TEST_CASE( picsar_equispaced_2d_table_setter)
{
    auto tab_2d = make_2d_table();
    const auto val = 1000.0;
    const int i = 10;
    const int j = 20;
    tab_2d.set_val(i,j,val);
    const auto xx = tab_2d.get_x_coord(i);
    const auto yy = tab_2d.get_y_coord(j);
    const auto res =tab_2d.interp(xx,yy);

    BOOST_CHECK_SMALL((res - val)/val, tolerance<double>());
}

// ***Test equispaced_2d_table set all
BOOST_AUTO_TEST_CASE( picsar_equispaced_2d_table_setall)
{
    auto tab_2d = make_2d_table();
    const auto val = 42.0;
    const auto vec = std::vector<double>(
        static_cast<std::size_t>(tab_2d.get_how_many_x())*
        static_cast<std::size_t>(tab_2d.get_how_many_y()),
        val);

    tab_2d.set_all_vals(vec);

    const auto x0 = tab_2d.get_x_min();
    const auto x1 = tab_2d.get_x_max();
    const auto x2 = 0.376*(x1-x0) + x0;

    const auto y0 = tab_2d.get_y_min();
    const auto y1 = tab_2d.get_y_max();
    const auto y2 = 0.688*(y1-y0) + y0;

    BOOST_CHECK_SMALL((tab_2d.interp(x0,y0) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x0,y1) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x0,y2) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x1,y0) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x1,y1) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x1,y2) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x2,y0) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x2,y1) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((tab_2d.interp(x2,y2) - val)/val, tolerance<double>());
}


// ***Test equispaced_2_table equality
BOOST_AUTO_TEST_CASE( picsar_equispaced_2d_table_equality)
{
    auto tab_2d = make_2d_table();
    auto tab_2d_2 = make_2d_table();
    BOOST_CHECK_EQUAL(tab_2d == tab_2d_2, true);

    tab_2d.set_val(0, 10.0);
    BOOST_CHECK_EQUAL(tab_2d == tab_2d_2, false);
}

// ***Test equispaced_2d_table serialization
BOOST_AUTO_TEST_CASE( picsar_equispaced_2d_table_serialization)
{
    auto tab_2d = make_2d_table();
    auto raw_data = tab_2d.serialize();
    auto tab_2d_2 = equispaced_2d_table<double, std::vector<double>>{raw_data};

    BOOST_CHECK_EQUAL(tab_2d.get_x_min(), tab_2d_2.get_x_min());
    BOOST_CHECK_EQUAL(tab_2d.get_x_max(), tab_2d_2.get_x_max());
    BOOST_CHECK_EQUAL(tab_2d.get_y_min(), tab_2d_2.get_y_min());
    BOOST_CHECK_EQUAL(tab_2d.get_y_max(), tab_2d_2.get_y_max());
    BOOST_CHECK_EQUAL(tab_2d.get_x_size(), tab_2d_2.get_x_size());
    BOOST_CHECK_EQUAL(tab_2d.get_y_size(), tab_2d_2.get_y_size());
    BOOST_CHECK_EQUAL(tab_2d.get_how_many_x(), tab_2d_2.get_how_many_x());
    BOOST_CHECK_EQUAL(tab_2d.get_how_many_y(), tab_2d_2.get_how_many_y());

    for(int i = 0; i < tab_2d.get_how_many_x()*tab_2d.get_how_many_x(); ++i){
        BOOST_CHECK_EQUAL(tab_2d.get_values_reference()[i],
            tab_2d_2.get_values_reference()[i]);
    }
}

// ***Test generic_2d_table constructor and getters
BOOST_AUTO_TEST_CASE( picsar_generic_2d_table_constructor_getters)
{
    auto gtab_2d = make_generic_2d_table();
    const auto const_gtab_2d = make_generic_2d_table();
    auto copy_gtab_2d = gtab_2d;

    check_generic_table_2d(gtab_2d);
    check_generic_table_2d(const_gtab_2d);
    check_generic_table_2d(copy_gtab_2d);
}

// ***Test generic_2d_table single coordinate interpolation
BOOST_AUTO_TEST_CASE( picsar_generic_2d_table_interp_one_coord)
{
    const auto const_gtab_2d = make_generic_2d_table();

    check_generic_table_2d_interp_one_coord(const_gtab_2d);
}

// ***Test generic_2d_table interp
BOOST_AUTO_TEST_CASE( picsar_generic_2d_table_constructor_interp)
{
    check_generic_table_2d_interp(0, 0, 1000.0);
    check_generic_table_2d_interp(gxsize-1, gysize-1, -1000.0);
    check_generic_table_2d_interp(gxsize-1, 0, 1000.0);
    check_generic_table_2d_interp(0, gysize-1, -1000.0);
    check_generic_table_2d_interp(5, 5, 1000.0);
    check_generic_table_2d_interp(gxsize-5, gysize-5, -1000.0);
    check_generic_table_2d_interp(gxsize-5, 5, 1000.0);
    check_generic_table_2d_interp(5, gysize-5, -1000.0);
}

// ***Test generic_2d_table equality
BOOST_AUTO_TEST_CASE( picsar_generic_2d_table_equality)
{
    auto gtab_2d = make_generic_2d_table();
    auto gtab_2d_2 = make_generic_2d_table();
    BOOST_CHECK_EQUAL(gtab_2d == gtab_2d_2, true);

    gtab_2d.set_val(0, 10.0);
    gtab_2d.set_val(gxsize/4, -7.0);
    gtab_2d.set_val(gxsize/2, 7.0);
    gtab_2d.set_val(gxsize-1, 100.0);
    BOOST_CHECK_EQUAL(gtab_2d == gtab_2d_2, false);
}

// ***Test generic_2d_table set all
BOOST_AUTO_TEST_CASE( picsar_generic_2d_table_setall)
{
    auto gtab_2d = make_generic_2d_table();
    const auto val = 42.0;
    const auto vec = std::vector<double>(
        static_cast<std::size_t>(tab_2d.get_how_many_x())*
        static_cast<std::size_t>(tab_2d.get_how_many_y()),
        val);

    gtab_2d.set_all_vals(vec);

    const auto x0 = tab_2d.get_x_min();
    const auto x1 = tab_2d.get_x_max();
    const auto x2 = 0.376*(x1-x0) + x0;

    const auto y0 = tab_2d.get_y_min();
    const auto y1 = tab_2d.get_y_max();
    const auto y2 = 0.688*(y1-y0) + y0;

    BOOST_CHECK_SMALL((gtab_2d.interp(x0,y0) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x0,y1) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x0,y2) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x1,y0) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x1,y1) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x1,y2) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x2,y0) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x2,y1) - val)/val, tolerance<double>());
    BOOST_CHECK_SMALL((gtab_2d.interp(x2,y2) - val)/val, tolerance<double>());
}

// ***Test generic_2d_table serialization
BOOST_AUTO_TEST_CASE( picsar_generic_2d_table_serialization)
{
    auto gtab_2d = make_generic_2d_table();
    auto raw_data = gtab_2d.serialize();
    auto gtab_2d_2 = Generic2DTableType{raw_data};

    BOOST_CHECK_EQUAL(gtab_2d.get_x_min(), gtab_2d_2.get_x_min());
    BOOST_CHECK_EQUAL(gtab_2d.get_x_max(), gtab_2d_2.get_x_max());
    BOOST_CHECK_EQUAL(gtab_2d.get_y_min(), gtab_2d_2.get_y_min());
    BOOST_CHECK_EQUAL(gtab_2d.get_y_max(), gtab_2d_2.get_y_max());
    BOOST_CHECK_EQUAL(gtab_2d.get_x_size(), gtab_2d_2.get_x_size());
    BOOST_CHECK_EQUAL(gtab_2d.get_y_size(), gtab_2d_2.get_y_size());
    BOOST_CHECK_EQUAL(gtab_2d.get_how_many_x(), gtab_2d_2.get_how_many_x());
    BOOST_CHECK_EQUAL(gtab_2d.get_how_many_y(), gtab_2d_2.get_how_many_y());

    for(int i = 0; i < gtab_2d.get_how_many_x()*gtab_2d.get_how_many_x(); ++i){
        BOOST_CHECK_EQUAL(gtab_2d.get_values_reference()[i],
            gtab_2d_2.get_values_reference()[i]);
    }

    const auto xt = 0.5*(gxmax - gxmin);
    const auto yt = 0.5*(gymax - gymin);

    BOOST_CHECK_EQUAL(gtab_2d.interp(xt, yt), gtab_2d_2.interp(xt, yt));
}
