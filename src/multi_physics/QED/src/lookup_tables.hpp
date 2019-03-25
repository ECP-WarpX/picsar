#ifndef __PICSAR_MULTIPHYSICS_LOOKUP_TABLES__
#define __PICSAR_MULTIPHYSICS_LOOKUP_TABLES__

//This .hpp file contais the implementation of 1D and 2D lookup tables.

#include<functional>
#include<array>
#include<vector>
#include<utility>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        //1D lookup table
        template<typename _REAL>
        class lookup_1d
        {
            public:
                //Default empty constructor
                lookup_1d () = default;

                //Constructor: requires coordinates and data
                //coordinates should be sorted!
                lookup_1d (std::vector<_REAL> coords, std::vector<_REAL> data);

                //Copy constructor
                lookup_1d(lookup_1d& other);

                //Move constructor
                lookup_1d(lookup_1d&& other);

                //Assignment operator
                lookup_1d&  operator= (const lookup_1d& );

                //Get a copy of the coordinates
                std::vector<_REAL> get_coords() const;

                //Get a reference to the coordinates
                std::vector<_REAL>& ref_coords();

                //Get a reference to the data
                std::vector<_REAL>& ref_data();

                //Check if the table is initialized
                bool is_init() const;

                //Linear equispaced interpolation
                PXRMP_FORCE_INLINE
                _REAL interp_linear_equispaced(_REAL where) const;

                //Linear interpolation
                PXRMP_FORCE_INLINE
                _REAL interp_linear(_REAL where) const;

            private:
                std::vector<_REAL> coords;
                std::vector<_REAL> data;

                bool init_flag = false;
        };

        //2D lookup table
        template<typename _REAL>
        class lookup_2d
        {
            public:
                //Default empty constructor
                lookup_2d () = default;

                //Constructor: requires coordinates and data
                //coordinates should be sorted!
                //Row major convention is used.
                lookup_2d(
                    std::array<std::vector<_REAL>,2> coords,
                    std::vector<_REAL> data);

                //Copy constructor
                lookup_2d(lookup_2d& other);

                //Move constructor
                lookup_2d(lookup_2d&& other);

                //Assignment operator
                lookup_2d&  operator= (const lookup_2d& );

                //Get a copy of the coordinates
                std::array<std::vector<_REAL>,2> get_coords();

                //Check if the table is initialized
                bool is_init() const;

                //Performs linear interpolation
                PXRMP_FORCE_INLINE
                _REAL interp_linear(_REAL where_x, _REAL where_y);


            private:
                std::array<std::vector<_REAL>,2> coords;
                //Underlying data is stored in a vector. Row major convention
                //is used.
                std::vector<_REAL> data;

                bool init_flag = false;

                //Row major accessor: converts coordinates indices to
                //the index
                PXRMP_FORCE_INLINE
                size_t row_major(
                size_t i, size_t j, size_t c1_size, size_t c2_size) const;
        };
    }
}

//############################################### Implementation

//_______________________1D table_______________________________

//Constructor: requires coordinates and data
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>::
lookup_1d
(std::vector<_REAL> coords, std::vector<_REAL> data):
    coords{coords}, data{data}
{
    init_flag = true;
}

//Copy constructor
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>::
lookup_1d(lookup_1d& other):
    coords{other.coords}, data{other.data}, init_flag{other.init_flag}
{}

//Move constructor
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>::
lookup_1d(lookup_1d&& other):
    coords{std::move(other.coords)}, data{std::move(other.data)},
    init_flag{other.init_flag}
{}

//Assignment operator
template<typename _REAL>
picsar::multi_physics::lookup_1d<_REAL>&
picsar::multi_physics::lookup_1d<_REAL>::
 operator= (const lookup_1d& other)
{
    if (this != &other){
        this->coords = other.coords;
        this->data = other.data;
        this->init_flag = other.init_flag;
    }
    return *this;
}

//Get a copy of the coordinates
template<typename _REAL>
std::vector<_REAL>
picsar::multi_physics::lookup_1d<_REAL>::
get_coords() const
{
    return coords;
}

template<typename _REAL>
std::vector<_REAL>&
picsar::multi_physics::lookup_1d<_REAL>::
ref_coords()
{
    return coords;
}

//Get a reference to the data
template<typename _REAL>
std::vector<_REAL>&
picsar::multi_physics::lookup_1d<_REAL>::
ref_data()
{
    return data;
}

//Checks if the table is initialized
template<typename _REAL>
bool
picsar::multi_physics::lookup_1d<_REAL>::
is_init() const
{
    return init_flag;
}

//Performs linear equispaced interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_1d<_REAL>::
interp_linear_equispaced(_REAL where) const
{
    _REAL xmin = coords.front();
    _REAL xmax = coords.back();
    _REAL yleft = data.front();
    _REAL yright = data.back();

    return yleft + ((where-xmin)/(xmax-xmin))*(yright-yleft);
}

//Performs linear interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_1d<_REAL>::
interp_linear(_REAL where) const
{
    auto iter_right = std::upper_bound(coords.begin(), coords.end(),where);
    size_t idx_right = std::distance(coords.begin(), iter_right);
    size_t idx_left = idx_right-1;

    _REAL xleft = coords[idx_left];
    _REAL xright = coords[idx_right];
    _REAL yleft = data[idx_left];
    _REAL yright = data[idx_right];

    return yleft + ((yright-yleft)/(xright-xleft))*(where-xleft);
}

//_______________________2D table_______________________________

//Constructor: requires coordinates and data
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>::
lookup_2d(
std::array<std::vector<_REAL>,2> coords,
std::vector<_REAL> data):
    coords{coords}, data{data}
{
    init_flag = true;
}

//Copy constructor
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>::
lookup_2d(lookup_2d& other):
    coords{other.coords}, data{other.data}, init_flag{other.init_flag}
{}

//Move constructor
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>::
lookup_2d(lookup_2d&& other):
    coords{std::move(other.coords)}, data{std::move(other.data)},
     init_flag{std::move(other.init_flag)}
{}

//Assignment operator
template<typename _REAL>
picsar::multi_physics::lookup_2d<_REAL>&
picsar::multi_physics::lookup_2d<_REAL>::
 operator= (const lookup_2d& other)
{
    if (this != &other){
        this->coords = other.coords;
        this->data = other.data;
        this->init_flag = other.init_flag;
    }
    return *this;
}



//Get a copy of the coordinates
template<typename _REAL>
std::array<std::vector<_REAL>,2>
picsar::multi_physics::lookup_2d<_REAL>::
get_coords()
{
    return coords;
}

//Checks if the table is initialized
template<typename _REAL>
bool
picsar::multi_physics::lookup_2d<_REAL>::
is_init() const
{
    return init_flag;
}

//Performs a linear interpolation
template<typename _REAL>
PXRMP_FORCE_INLINE
_REAL
picsar::multi_physics::lookup_2d<_REAL>::
interp_linear(_REAL where_x, _REAL where_y)
{
    auto it_x_right =
        std::upper_bound(coords[0].begin(), coords[0].end(),where_x);
    auto it_y_right =
        std::upper_bound(coords[1].begin(), coords[1].end(),where_y);
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

    _REAL zll = data[row_major(idx_x_left,idx_y_left,sx,sy)];
    _REAL zlr = data[row_major(idx_x_left,idx_y_right,sx,sy)];
    _REAL zrl = data[row_major(idx_x_right,idx_y_left,sx,sy)];
    _REAL zrr = data[row_major(idx_x_right,idx_y_right,sx,sy)];

    _REAL wll = (xright - where_x)*(yright -where_y);
    _REAL wlr = (xright - where_x)*(where_y - yleft);
    _REAL wrl = (where_x - xleft)*(yright -where_y);
    _REAL wrr = (where_x - xleft)*(where_y - yleft);

    _REAL w_norm = (xright-xleft)*(yright-yleft);

    return (wll*zll + wlr*zlr + wrl*zrl + wrr*zrr)/w_norm;
}

//Row major access to underlying data
template<typename _REAL>
PXRMP_FORCE_INLINE
size_t
picsar::multi_physics::lookup_2d<_REAL>::
row_major(
size_t i, size_t j, size_t, size_t c2_size) const
{
    return i*c2_size + j;
}
//______________________
#endif// __PICSAR_MULTIPHYSICS_LOOKUP_TABLES__
