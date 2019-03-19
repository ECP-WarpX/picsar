#ifndef __PICSAR_MULTIPHYSICS_LOOKUP_TABLES__
#define __PICSAR_MULTIPHYSICS_LOOKUP_TABLES__

//This .hpp file contais the implementation of 1D and 2D lookup tables.
//The interpolator function should be provided when the lookup table
//is initilialized.

#include<functional>
#include<array>
#include<vector>
#include<utility>

//Should be included by all the src files of the library
#include "qed_commons.h"


//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //Some useful aliases
        template<typename _REAL>
        using interpolator_1d =
        std::function<_REAL(_REAL val,
        const std::vector<_REAL>& coords,
        const std::vector<_REAL>& data)>;

        using accessor_2d =
        std::function<size_t(size_t i, size_t j, size_t c1_size, size_t c2_size)>;

        template<typename _REAL>
        using interpolator_2d =
        std::function<_REAL(_REAL valx, _REAL valy,
        const std::array<std::vector<_REAL>,2>& coords,
        const std::vector<_REAL>& data, accessor_2d accessor)>;

        //1D lookup table
        template<typename _REAL>
        class lookup_1d
        {
            public:
                //Constructor: requires coordinates, data and an interpolator
                //interpolator should be a lambda
                //f(_REAL val, const vector<_REAL>& coords, const vector<_REAL>& data)-> _REAL .
                //coordinates should be sorted!
                lookup_1d (std::vector<_REAL> coords,
                std::vector<_REAL> data, interpolator_1d<_REAL> inteporlator);

                //Copy constructor
                lookup_1d(lookup_1d& other);

                //Move constructor
                lookup_1d(lookup_1d&& other);

                //Get a copy of the coordinates
                std::vector<_REAL> get_coords();

                //Performs the interpolation
                PXRMP_FORCE_INLINE
                _REAL interp(_REAL where);

                //----Static pre-defined interpolator functions

                //Linear equispaced interpolation
                PXRMP_FORCE_INLINE
                static _REAL linear_equispaced_interpolation(
                _REAL val,
                const std::vector<_REAL>& coords,
                const std::vector<_REAL>& data);

                //Linear interpolation
                PXRMP_FORCE_INLINE
                static _REAL linear_interpolation(
                _REAL val,
                const std::vector<_REAL>& coords,
                const std::vector<_REAL>& data);
                //______________________

            private:
                std::vector<_REAL> coords;
                std::vector<_REAL> data;

                interpolator_1d<_REAL> interpolator;
        };

        //2D lookup table
        template<typename _REAL>
        class lookup_2d
        {
            public:
                //Constructor: requires coordinates, data, an interpolator and an accessor
                //interpolator should be a lambda
                //f(array<REAL,2>, array<vector<_REAL>,2>& coords, vector<_REAL>& data)-> _REAL .
                //coordinates should be sorted!
                lookup_2d(
                    std::array<std::vector<_REAL>,2> coords,
                    std::vector<_REAL> data,
                    interpolator_2d<_REAL> inteporlator,
                    accessor_2d accessor);

                //Copy constructor
                lookup_2d(lookup_2d& other);

                //Move constructor
                lookup_2d(lookup_2d&& other);

                //Get a copy of the coordinates
                std::array<std::vector<_REAL>,2> get_coords();

                //Performs the interpolation
                PXRMP_FORCE_INLINE
                _REAL interp(_REAL where_x, _REAL where_y);

                //----Static pre-defined interpolator&accessor functions

                //Linear interpolation
                PXRMP_FORCE_INLINE
                static _REAL linear_interpolation(
                _REAL valx, _REAL valy,
                const std::array<std::vector<_REAL>,2>& coords,
                const std::vector<_REAL>& data, accessor_2d accessor);

                //Row major accessor
                PXRMP_FORCE_INLINE
                static size_t row_major(
                size_t i, size_t j, size_t c1_size, size_t);
                //______________________


            private:
                std::array<std::vector<_REAL>,2> coords;
                std::vector<_REAL> data;

                interpolator_2d<_REAL>  interpolator;
                accessor_2d accessor;
        };
    }
}

//############################################### Implementation

//_______________________1D table_______________________________

//Constructor: requires coordinates, data and an interpolator
//interpolator should be a lambda
//f(array<REAL,2>, array<vector<_REAL>,2>& coords, vector<_REAL>& data)-> _REAL .
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>::
lookup_1d
(std::vector<_REAL> coords, std::vector<_REAL> data,
interpolator_1d<_REAL> interpolator):
    coords{coords}, data{data}, interpolator{interpolator}
{}

//Copy constructor
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>::
lookup_1d(lookup_1d& other):
    coords{other.coords}, data{other.data}, interpolator{other.interpolator}
{}

//Move constructor
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>::
lookup_1d(lookup_1d&& other):
    coords{std::move(other.coords)}, data{std::move(other.data)},
    interpolator{std::move(other.interpolator)}
{}

//Get a copy of the coordinates
template<typename _REAL>
std::vector<_REAL>
picsar::multi_physics::lookup_1d<_REAL>::
get_coords()
{
    return coords;
}

//Performs the interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_1d<_REAL>::
interp(_REAL where)
{
    return interpolator(where, coords, data);
}

//----Static pre-defined interpolator functions
//Linear equispaced interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_1d<_REAL>::
linear_equispaced_interpolation(
_REAL val,
const std::vector<_REAL>& coords,
const std::vector<_REAL>& data)
{
    _REAL xmin = coords.front();
    _REAL xmax = coords.back();
    _REAL yleft = data.front();
    _REAL yright = data.back();

    return yleft + ((val-xmin)/(xmax-xmin))*(yright-yleft);
}


//Linear interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_1d<_REAL>::
linear_interpolation(
_REAL val,
const std::vector<_REAL>& coords,
const std::vector<_REAL>& data)
{
    auto iter_right = std::upper_bound(coords.begin(), coords.end(),val);
    size_t idx_right = std::distance(coords.begin(), iter_right);
    size_t idx_left = idx_right-1;

    _REAL xleft = coords[idx_left];
    _REAL xright = coords[idx_right];
    _REAL yleft = data[idx_left];
    _REAL yright = data[idx_right];

    return yleft + ((yright-yleft)/(xright-xleft))*(val-xleft);
}
//_______________________2D table_______________________________

//Constructor: requires coordinates, data and an interpolator
//interpolator should be a lambda
//f(array<REAL,2>, array<vector<_REAL>,2>& coords, vector<_REAL>& data)-> _REAL.
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>::
lookup_2d(
std::array<std::vector<_REAL>,2> coords,
std::vector<_REAL> data,
interpolator_2d<_REAL> interpolator,
accessor_2d accessor):
    coords{coords}, data{data}, interpolator{interpolator}, accessor{accessor}
    {}

//Copy constructor
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>::
lookup_2d(lookup_2d& other):
    coords{other.coords}, data{other.data}, interpolator{other.interpolator},
    accessor{other.accessor}
{}

//Move constructor
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>::
lookup_2d(lookup_2d&& other):
    coords{std::move(other.coords)}, data{std::move(other.data)},
    interpolator{std::move(other.interpolator)},
    accessor{std::move(other.accessor)}
{}


//Get a copy of the coordinates
template<typename _REAL>
std::array<std::vector<_REAL>,2>
picsar::multi_physics::lookup_2d<_REAL>::
get_coords()
{
    return coords;
}

//Performs the interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_2d<_REAL>::
interp(_REAL where_x, _REAL where_y)
{
    return interpolator(where_x, where_y, coords, data, accessor);
}

//----Static pre-defined interpolator functions
//Linear interpolation in 2D
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_2d<_REAL>::
linear_interpolation(
_REAL valx, _REAL valy,
const std::array<std::vector<_REAL>,2>& coords,
const std::vector<_REAL>& data,
picsar::multi_physics::accessor_2d accessor)
{
    auto it_x_right = std::upper_bound(coords[0].begin(), coords[0].end(),valx);
    auto it_y_right = std::upper_bound(coords[1].begin(), coords[1].end(),valy);
    size_t idx_x_right = std::distance(coords[0].begin(), it_x_right);
    size_t idx_y_right = std::distance(coords[1].begin(), it_y_right);
    size_t idx_x_left = idx_x_right - 1;
    size_t idx_y_left = idx_y_right - 1;

    _REAL xleft = coords[0][idx_x_left];
    _REAL xright = coords[0][idx_x_right];
    _REAL yleft = coords[1][idx_y_left];
    _REAL yright = coords[1][idx_y_right];

    size_t sx = coords[0].size();
    size_t sy = coords[1].size();

    _REAL zll = data[accessor(idx_x_left,idx_y_left,sx,sy)];
    _REAL zlr = data[accessor(idx_x_left,idx_y_right,sx,sy)];
    _REAL zrl = data[accessor(idx_x_right,idx_y_left,sx,sy)];
    _REAL zrr = data[accessor(idx_x_right,idx_y_right,sx,sy)];

    _REAL wll = (xright - valx)*(yright -valy);
    _REAL wlr = (xright - valx)*(valy - yleft);
    _REAL wrl = (valx - xleft)*(yright -valy);
    _REAL wrr = (valx - xleft)*(valy - yleft);

    _REAL w_norm = (xright-xleft)*(yright-yleft);

    return (wll*zll + wlr*zlr + wrl*zrl + wrr*zrr)/w_norm;
}

//Row major access to underlying data
template<typename _REAL>
PXRMP_FORCE_INLINE
size_t
picsar::multi_physics::lookup_2d<_REAL>::
row_major(
size_t i, size_t j, size_t c1_size, size_t )
{
    return i*c1_size + j;
}
//______________________
#endif// __PICSAR_MULTIPHYSICS_LOOKUP_TABLES__
