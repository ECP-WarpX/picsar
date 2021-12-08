#ifndef PICSAR_TABLES
#define PICSAR_TABLES

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

// Interpolation is used
#include "picsar_qed/utils/picsar_algo.hpp"
// Serialization is used to export tables to byte vectors
// and to read tables from byte arrays
#include "picsar_qed/utils/serialization.hpp"
// Some mathematical constants are used
#include "picsar_qed/math/math_constants.h"
//Uses floor
#include "picsar_qed/math/cmath_overloads.hpp"

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace picsar{
namespace multi_physics{
namespace containers{

    //________________ 1D equispaced table _____________________________________

    /**
    * This class implements a generic equispaced 1D lookup table
    * for a function f(x)
    *
    * @tparam RealType the floating point type to be used (e.g. double or float)
    * @tparam VectorType the vector type to be used (e.g. std::vector<double>)
    */
    template <typename RealType, class VectorType>
    class equispaced_1d_table{

    public:

        /**
        * Constructor (not usable on GPUs)
        *
        * @param[in] x_min the leftmost extreme of the x coordinates
        * @param[in] x_max the rightmost extreme of the x coordinates
        * @param[in] values the values of the function that should be interpolated.
        */
        equispaced_1d_table(
            RealType x_min, RealType x_max, VectorType values):
            m_x_min{x_min}, m_x_max{x_max}, m_values{values}
            {
                m_how_many_x = static_cast<int>(values.size());
                m_x_size = x_max - x_min;
                m_dx = m_x_size/(m_how_many_x-1);
            }

        /**
        * Empty constructor
        */
        constexpr
        equispaced_1d_table(){}

        /**
        * Constructor from byte array (not usable on GPUs)
        *
        * @param[in] raw_data a const reference to a byte vector
        */
        equispaced_1d_table(const std::vector<char>& raw_data)
        {
            using namespace picsar::multi_physics::utils;

            constexpr size_t min_size =
                sizeof(char)+//double or float
                sizeof(m_x_min)+sizeof(m_x_max)+
                sizeof(m_how_many_x) + sizeof(m_dx);

            if (raw_data.size() < min_size)
                throw std::runtime_error("Binary data is too small \
                to be a 1D table.");

            auto it_raw_data = raw_data.begin();

            if (serialization::get_out<char>(it_raw_data) !=
                static_cast<char>(sizeof(RealType))){
                throw std::runtime_error("Mismatch between RealType \
                used to write and to read the 1D table");
            }

            m_x_min = serialization::get_out<decltype(m_x_min)>(it_raw_data);
            m_x_max = serialization::get_out<decltype(m_x_max)>(it_raw_data);
            m_x_size = m_x_max - m_x_min;
            if(m_x_size < 0)
                throw std::runtime_error("raw_data contains invalid data.");

            m_how_many_x =
                serialization::get_out<decltype(m_how_many_x)>(it_raw_data);
            m_dx = serialization::get_out<decltype(m_dx)>(it_raw_data);
            if(m_how_many_x <= 0)
                throw std::runtime_error("raw_data contains invalid data.");
            m_values = VectorType(m_how_many_x);
            auto vals =
                serialization::get_n_out<RealType>(it_raw_data, m_how_many_x);
            std::copy(vals.begin(), vals.end(), m_values.begin());
        }

        /**
        * Operator ==
        *
        * @param[in] rhs a const reference to a 1D table of the same type
        * @return true if rhs and *this are equal. false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator== (
            const equispaced_1d_table<RealType, VectorType> &rhs) const
        {
            return
                (m_x_min == rhs.m_x_min) &&
                (m_x_max == rhs.m_x_max) &&
                (m_x_size == rhs.m_x_size) &&
                (m_how_many_x == rhs.m_how_many_x) &&
                (m_dx == rhs.m_dx) &&
                (m_values == rhs.m_values);
        }

        /**
        * Returns the i-th coodinate along x axis
        *
        * @param[in] i the index of the desired coordinate
        * @return the i-th coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_coord(const int i) const noexcept
        {
                return i*m_dx + m_x_min;
        }

        /**
        * Returns the i-th value
        *
        * @param[in] i the index of the desired value
        * @return the i-th value
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_val(const size_t i) const noexcept
        {
                return m_values[i];
        }

        /**
        * Returns the number of points along x
        *
        * @return the number of points along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int get_how_many_x() const noexcept
        {
            return m_how_many_x;
        }

        /**
        * Returns the minimum x coordinate
        *
        * @return the minimum x coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_min() const noexcept
        {
            return m_x_min;
        }

        /**
        * Returns the maximum x coordinate
        *
        * @return the maximum x coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_max() const noexcept
        {
            return m_x_max;
        }

        /**
        * Returns the size along x
        *
        * @return the size along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_size() const noexcept
        {
            return m_x_size;
        }

        /**
        * Returns the size of the steps along x
        *
        * @return the size of the steps along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_dx() const noexcept
        {
            return m_dx;
        }

        /**
        * Returns an std::vector containing all the coordinates (not usable on
        * GPUs).
        *
        * @return all the coordinates
        */
        std::vector<RealType> get_all_coordinates()  const noexcept
        {
            auto all_coords = std::vector<RealType>(m_how_many_x);
            for (int i = 0; i < m_how_many_x; ++i){
                all_coords[i] = get_x_coord(i);
            }
            return all_coords;
        }

        /**
        * Returns a const reference to the underlying Vector holding value data
        *
        * @return a const reference to the underlying Vector holding value data
        */
        const VectorType& get_values_reference() const noexcept
        {
            return m_values;
        }

        /**
        * Perfoms a linear interpolation at a given position
        * Warning: where_x must be within [x_min, x_max]. No safety check
        * is performed!
        *
        * @param[in] where_x where to perform the interpolation
        * @return the interpolated value
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp(const RealType where_x) const noexcept
        {
            using namespace picsar::multi_physics::math;

            const auto idx_left = static_cast<int>(
                m_floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            if (idx_left == (m_how_many_x-1))
                return  m_values[m_how_many_x-1];
            const auto idx_right = idx_left + 1;

            const auto xleft = get_x_coord(idx_left);
            const auto xright = get_x_coord(idx_right);
            const auto yleft = m_values[idx_left];
            const auto yright = m_values[idx_right];

            return utils::linear_interp(xleft, xright, yleft, yright, where_x);
        }

        /**
        * Overwrites the i-th value (not designed for GPU usage)
        * Warning: no safety check on i is performed!
        *
        * @param[in] i index of the value to be overwritten
        * @param[in] what value to be written at position i
        */
        PXRMP_FORCE_INLINE
        void set_val(int i, RealType what)
        {
            m_values[i] = what;
        }

        /**
        * Returns a byte vector containing all the raw data of the
        * table (not usable on GPUs)
        *
        * @return a byte vector containing table data
        */
        std::vector<char> serialize() const
        {
            auto raw_data = std::vector<char>{};

            utils::serialization::put_in(
                static_cast<char>(sizeof(RealType)), raw_data);
            utils::serialization::put_in(m_x_min, raw_data);
            utils::serialization::put_in(m_x_max, raw_data);
            utils::serialization::put_in(m_how_many_x, raw_data);
            utils::serialization::put_in(m_dx, raw_data);
            for (auto val : m_values)
                utils::serialization::put_in(val, raw_data);

            return raw_data;
        }

    protected:
        RealType m_x_min = static_cast<RealType>(0.0); /* minimum x coordinate */
        RealType m_x_max = static_cast<RealType>(0.0);/* maximum x coordinate */
        RealType m_x_size = static_cast<RealType>(0.0); /* size along x */
        int m_how_many_x = 0; /* how many values are stored in the table */
        RealType m_dx = static_cast<RealType>(0.0); /* size of the step along x */
        VectorType m_values; /* values f(x) */
    };

    //__________________________________________________________________________

    //________________ 2D equispaced table _____________________________________

    /**
    * This class implements a generic equispaced 2D lookup table
    * for a function f(x, y)
    *
    * @tparam RealType the floating point type to be used (e.g. double or float)
    * @tparam VectorType the vector type to be used (e.g. std::vector<double>)
    */
    template <typename RealType, class VectorType>
    class equispaced_2d_table
    {
    public:

        /**
        * Constructor (not usable on GPUs)
        *
        * @param[in] x_min the leftmost extreme of the x coordinates
        * @param[in] x_max the rightmost extreme of the x coordinates
        * @param[in] y_min the leftmost extreme of the y coordinates
        * @param[in] y_max the rightmost extreme of the y coordinates
        * @param[in] how_many_x number of grid points along x
        * @param[in] how_many_y number of grid points along y
        * @param[in] values the values of the function that should be interpolated (row major order)
        */
        equispaced_2d_table(
            RealType x_min, RealType x_max,
            RealType y_min, RealType y_max,
            int how_many_x, int how_many_y,
            VectorType values):
            m_x_min{x_min}, m_x_max{x_max},
            m_y_min{y_min}, m_y_max{y_max},
            m_how_many_x{how_many_x}, m_how_many_y{how_many_y},
            m_values{values}
            {
                m_x_size = x_max - x_min;
                m_y_size = y_max - y_min;
                m_dx = m_x_size/(m_how_many_x-1);
                m_dy = m_y_size/(m_how_many_y-1);
            }

        /**
        * Empty constructor
        */
        constexpr
        equispaced_2d_table(){}

        /**
        * Constructor from byte array (not usable on GPUs)
        *
        * @param[in] raw_data a const reference to a byte vector
        */
        equispaced_2d_table(const std::vector<char>& raw_data)
        {
            using namespace picsar::multi_physics::utils;

            constexpr size_t min_size =
                sizeof(char)+//double or float
                sizeof(m_x_min)+sizeof(m_x_max)+
                sizeof(m_y_min)+sizeof(m_y_max)+
                sizeof(m_how_many_x)+
                sizeof(m_how_many_y)+
                sizeof(m_dx) + sizeof(m_dy);

            if (raw_data.size() < min_size)
                throw std::runtime_error("Binary data is too small \
                to be a 2D table.");

            auto it_raw_data = raw_data.begin();

            if (serialization::get_out<char>(it_raw_data) !=
                static_cast<char>(sizeof(RealType))){
                throw std::runtime_error("Mismatch between RealType \
                used to write and to read the 1D table");
            }

            m_x_min = serialization::get_out<decltype(m_x_min)>(it_raw_data);
            m_x_max = serialization::get_out<decltype(m_x_max)>(it_raw_data);
            m_y_min = serialization::get_out<decltype(m_y_min)>(it_raw_data);
            m_y_max = serialization::get_out<decltype(m_y_max)>(it_raw_data);
            m_x_size = m_x_max - m_x_min;
            m_y_size = m_y_max - m_y_min;
            if(m_x_size < 0)
                throw std::runtime_error("raw_data contains invalid data.");
            if(m_y_size < 0)
                throw std::runtime_error("raw_data contains invalid data.");

            m_how_many_x =
                serialization::get_out<decltype(m_how_many_x)>(it_raw_data);
            m_how_many_y =
                serialization::get_out<decltype(m_how_many_y)>(it_raw_data);
            m_dx = serialization::get_out<decltype(m_dx)>(it_raw_data);
            m_dy = serialization::get_out<decltype(m_dy)>(it_raw_data);
            if(m_how_many_x <= 0)
                throw std::runtime_error("raw_data contains invalid data.");
            if(m_how_many_y <= 0)
                throw std::runtime_error("raw_data contains invalid data.");
            m_values = VectorType(m_how_many_x*m_how_many_y);
            auto vals = serialization::get_n_out<RealType>(
                    it_raw_data,
                    m_how_many_x*m_how_many_y);
            std::copy(vals.begin(), vals.end(), m_values.begin());
        }

        /**
        * Operator ==
        *
        * @param[in] rhs a const reference to a 2D table of the same type
        * @return true if rhs and *this are equal. false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator== (
            const equispaced_2d_table<RealType, VectorType> &rhs) const
        {
            return
                (m_x_min == rhs.m_x_min) &&
                (m_x_max == rhs.m_x_max) &&
                (m_y_min == rhs.m_y_min) &&
                (m_y_max == rhs.m_y_max) &&
                (m_x_size == rhs.m_x_size) &&
                (m_y_size == rhs.m_y_size) &&
                (m_how_many_x == rhs.m_how_many_x) &&
                (m_how_many_y == rhs.m_how_many_y) &&
                (m_dx == rhs.m_dx) &&
                (m_dy == rhs.m_dy) &&
                (m_values == rhs.m_values);
        }

        /**
        * Returns the i-th coodinate along x axis
        *
        * @param[in] i the index of the desired coordinate
        * @return the i-th coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_coord(int i) const noexcept
        {
            return i*m_dx + m_x_min;
        }

        /**
        * Returns the j-th coodinate along y axis
        *
        * @param[in] j the index of the desired coordinate
        * @return the j-th coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_coord(int j) const noexcept
        {
            return j*m_dy + m_y_min;
        }

        /**
        * Returns the function value at (i, j)
        *
        * @param[in] i the index of the desired value along x
        * @param[in] j the index of the desired value along y
        * @return the value at (i,j)
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_val(int i, int j) const noexcept
        {
            return m_values[idx(i, j)];
        }

        /**
        * Returns the number of points along x
        *
        * @return the number of points along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int get_how_many_x() const noexcept
        {
            return m_how_many_x;
        }

        /**
        * Returns the minimum x coordinate
        *
        * @return the minimum x coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_min() const noexcept
        {
            return m_x_min;
        }

        /**
        * Returns the maximum x coordinate
        *
        * @return the maximum x coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_max() const noexcept
        {
            return m_x_max;
        }

        /**
        * Returns the size along x
        *
        * @return the size along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_size() const noexcept
        {
            return m_x_size;
        }

        /**
        * Returns the size of the steps along x
        *
        * @return the size of the steps along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_dx() const noexcept
        {
            return m_dx;
        }

        /**
        * Returns the number of points along y
        *
        * @return the number of points along y
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int get_how_many_y() const noexcept
        {
            return m_how_many_y;
        }

        /**
        * Returns the minimum y coordinate
        *
        * @return the minimum y coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_min() const noexcept
        {
            return m_y_min;
        }

        /**
        * Returns the maximum y coordinate
        *
        * @return the maximum y coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_max() const noexcept
        {
            return m_y_max;
        }

        /**
        * Returns the size along y
        *
        * @return the size along y
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_size() const noexcept
        {
            return m_y_size;
        }

        /**
        * Returns the size of the steps along y
        *
        * @return the size of the steps along y
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_dy() const noexcept
        {
            return m_dy;
        }

        /**
        * Returns an std::vector containing all the coordinates (not usable on
        * GPUs). Row major order is used.
        *
        * @return all the coordinates
        */
        std::vector<std::array<RealType,2>> get_all_coordinates() const noexcept
        {
            auto all_coords = std::vector<std::array<RealType,2>>(
                m_how_many_x*m_how_many_y, {0,0});
            int count = 0;
            for (int i = 0; i < m_how_many_x; ++i){
                for (int j = 0; j < m_how_many_y; ++j){
                    all_coords[count][0] = get_x_coord(i);
                    all_coords[count][1] = get_y_coord(j);
                    count++;
                }
            }
            return all_coords;
        }

        /**
        * Returns a const reference to the underlying Vector holding value data
        *
        * @return a const reference to the underlying Vector holding value data
        */
        const VectorType& get_values_reference() const noexcept
        {
            return m_values;
        }

        /**
        * Performs an interpolation a bilinear interpolation at given position.
        *
        * @param[in] where_x the position along x
        * @param[in] where_y the position along y
        * @return the result of the interpolation
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp(const RealType where_x, const RealType where_y) const noexcept
        {
            using namespace picsar::multi_physics::math;

            auto idx_x_left = static_cast<int>(
                m_floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            if (idx_x_left == (m_how_many_x-1))
                idx_x_left = m_how_many_x-2;
            const auto idx_x_right = idx_x_left + 1;
            auto idx_y_left = static_cast<int>(
                m_floor((m_how_many_y-1)*(where_y-m_y_min)/m_y_size));
            if (idx_y_left == (m_how_many_y-1))
                idx_y_left = m_how_many_y-2;
            const auto idx_y_right = idx_y_left + 1;

            const auto xleft = get_x_coord(idx_x_left);
            const auto xright = get_x_coord(idx_x_right);
            const auto yleft = get_y_coord(idx_y_left);
            const auto yright = get_y_coord(idx_y_right);
            const auto f_xl_yl = m_values[idx(idx_x_left, idx_y_left)];
            const auto f_xl_yr = m_values[idx(idx_x_left, idx_y_right)];
            const auto f_xr_yl = m_values[idx(idx_x_right, idx_y_left)];
            const auto f_xr_yr = m_values[idx(idx_x_right, idx_y_right)];

            return utils::bilinear_interp(
                xleft, xright, yleft, yright,
                f_xl_yl, f_xl_yr, f_xr_yl, f_xr_yr,
                where_x, where_y);
        }

        /**
        * Performs an interpolation a linear interpolation at given position,
        * in the special case in which the second coordinate is exactly a grid
        * point along y
        *
        * @param[in] where_x the position along x
        * @param[in] j the index of the position along y
        * @return the result of the interpolation
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp_first_coord(RealType where_x, int j) const noexcept
        {
            using namespace picsar::multi_physics::math;

            auto idx_left = static_cast<int>(
                m_floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            if (idx_left == (m_how_many_x-1))
                idx_left = m_how_many_x-2;
            const auto idx_right = idx_left + 1;

            const auto xleft = idx_left*m_dx + m_x_min;
            const auto xright = idx_right*m_dx + m_x_min;
            const auto left_val = m_values[idx(idx_left,j)];
            const auto right_val = m_values[idx(idx_right,j)];

            return utils::linear_interp(xleft, xright, left_val, right_val, where_x);
        }

        /**
        * Performs an interpolation a linear interpolation at given position,
        * in the special case in which the first coordinate is exactly a grid
        * point along x
        *
        * @param[in] i the position along x
        * @param[in] where_y the position along y
        * @return the result of the interpolation
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp_second_coord(int i, RealType where_y) const noexcept
        {
            using namespace picsar::multi_physics::math;

            auto idx_left = static_cast<int>(
                m_floor((m_how_many_y-1)*(where_y-m_y_min)/m_y_size));
            if (idx_left == (m_how_many_y-1))
                idx_left = m_how_many_y-2;
            const auto idx_right = idx_left + 1;

            const auto left_val = m_values[idx(i, idx_left)];
            const auto right_val = m_values[idx(i, idx_right)];
            const auto yleft = idx_left*m_dy + m_y_min;
            const auto yright = idx_right*m_dy + m_y_min;

            return utils::linear_interp(yleft, yright, left_val, right_val, where_y);
        }

        /**
        * Overwrites the (i,j) value (not designed for GPU usage)
        * Warning: no safety check on i and j is performed!
        *
        * @param[in] i index of the value to be overwritten
        * @param[in] j index of the value to be overwritten
        * @param[in] what value to be written at position (i,j)
        */
        PXRMP_FORCE_INLINE
        void set_val(const int i, const int j, const RealType what)
        {
            m_values[idx(i, j)] = what;
        }

        /**
        * Overwrites the i-th value (not designed for GPU usage)
        * Warning: no safety check on i is performed!
        *
        * @param[in] i index of the value to be overwritten
        * @param[in] what value to be written at position i
        */
        PXRMP_FORCE_INLINE
        void set_val(const int i, const RealType what)
        {
            m_values[i] = what;
        }

        /**
        * Returns a byte vector containing all the raw data of the
        * table (not usable on GPUs)
        *
        * @return a byte vector containing table data
        */
        std::vector<char> serialize() const
        {
            auto raw_data = std::vector<char>{};

            utils::serialization::put_in(
                static_cast<char>(sizeof(RealType)), raw_data);
            utils::serialization::put_in(m_x_min, raw_data);
            utils::serialization::put_in(m_x_max, raw_data);
            utils::serialization::put_in(m_y_min, raw_data);
            utils::serialization::put_in(m_y_max, raw_data);
            utils::serialization::put_in(m_how_many_x, raw_data);
            utils::serialization::put_in(m_how_many_y, raw_data);
            utils::serialization::put_in(m_dx, raw_data);
            utils::serialization::put_in(m_dy, raw_data);
            for (auto val : m_values)
                utils::serialization::put_in(val, raw_data);

            return raw_data;
        }

    protected:

            RealType m_x_min = static_cast<RealType>(0.0); /* minimum x coordinate */
            RealType m_x_max = static_cast<RealType>(0.0); /* maximum x coordinate */
            RealType m_y_min = static_cast<RealType>(0.0); /* minimum y coordinate */
            RealType m_y_max = static_cast<RealType>(0.0); /* maximum y coordinate */
            RealType m_x_size = static_cast<RealType>(0.0); /* size along x */
            RealType m_y_size = static_cast<RealType>(0.0); /* size along y */
            int m_how_many_x = 0; /* how many grid points along x */
            int m_how_many_y = 0; /* how many grid points along y */
            RealType m_dx = static_cast<RealType>(0.0); /* step size along x */
            RealType m_dy = static_cast<RealType>(0.0); /* step size along y */
            VectorType m_values; /* values f(x,y) */

        /**
        * Function values are stored internally in a 1D vector.
        * This function calculates the index along this vector
        * from the coordinate indices i,j. Row major order is used.
        *
        * @param[in] i index along x
        * @param[in] j index along y
        * @return index along internal 1D vector
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int idx(int i, int j) const noexcept
        {
            return i*m_how_many_y + j;
        }
    };

    //__________________________________________________________________________

    //________________ 2D generic table _____________________________________

    /**
    * This class implements a generic 2D lookup table for a function f(x, y),
    * using functors to map indices to coordinates and vice versa.
    * Functors must be trivially copiable and must implement a serialize() method
    * and a deserialize(CharIterator&) static method (see /physics/quantum_sync/quantum_sync_engine_tables_detail.hpp
    * for an example).
    *
    * @tparam RealType the floating point type to be used (e.g. double or float)
    * @tparam VectorType the vector type to be used (e.g. std::vector<double>)
    * @tparam XMapFunctor a functor to map the index on the x axis to the corresponding coordinate x
    * @tparam YMapFunctor a functor to map the index on the y axis to the corresponding coordinate y
    * @tparam IXMapFunctor a functor to map the coordinate x to the corresponding index
    * @tparam IYMapFunctor a functor to map the coordinate y to the corresponding index
    */
    template <
        typename RealType, typename VectorType,
        typename XMapFunctor, typename YMapFunctor,
        typename IXMapFunctor, typename IYMapFunctor>
    class generic_2d_table
    {
    public:

        /**
        * Constructor
        *
        * @param[in] how_many_x number of grid points along x
        * @param[in] how_many_y number of grid points along y
        * @param[in] values the values of the function that should be interpolated (row major order)
        * @param[in] x_map_functor the functor mapping the index on the x axis to the corresponding coordinate x
        * @param[in] y_map_functor the functor mapping the index on the y axis to the corresponding coordinate y
        * @param[in] ix_map_functor the functor mapping the coordinate x to the corresponding index
        * @param[in] iy_map_functor the functor mapping the coordinate y to the corresponding index
        */
        generic_2d_table(
            int how_many_x, int how_many_y,
            VectorType values,
            XMapFunctor x_map_functor,
            YMapFunctor y_map_functor,
            IXMapFunctor ix_map_functor,
            IYMapFunctor iy_map_functor):
            m_how_many_x{how_many_x}, m_how_many_y{how_many_y},
            m_values{values},
            m_x_map_functor{x_map_functor}, m_y_map_functor{y_map_functor},
            m_ix_map_functor{ix_map_functor}, m_iy_map_functor{iy_map_functor}
            {}

        /**
        * Empty constructor
        */
        constexpr
        generic_2d_table(){}

        /**
        * Constructor from byte array (not usable on GPUs)
        *
        * @param[in] raw_data a const reference to a byte vector
        */
        generic_2d_table(const std::vector<char>& raw_data)
        {
            using namespace picsar::multi_physics::utils;

            constexpr size_t min_size =
                sizeof(char)+//double or float
                sizeof(m_how_many_x)+
                sizeof(m_how_many_y);

            if (raw_data.size() < min_size)
                throw std::runtime_error("Binary data is too small \
                to be a 2D table.");

            auto it_raw_data = raw_data.begin();

            if (serialization::get_out<char>(it_raw_data) !=
                static_cast<char>(sizeof(RealType))){
                throw std::runtime_error("Mismatch between RealType \
                used to write and to read the 1D table");
            }

            m_how_many_x =
                serialization::get_out<decltype(m_how_many_x)>(it_raw_data);
            m_how_many_y =
                serialization::get_out<decltype(m_how_many_y)>(it_raw_data);
            if(m_how_many_x <= 0)
                throw std::runtime_error("raw_data contains invalid data.");
            if(m_how_many_y <= 0)
                throw std::runtime_error("raw_data contains invalid data.");
            m_values = VectorType(m_how_many_x*m_how_many_y);
            auto vals = serialization::get_n_out<RealType>(
                    it_raw_data,
                    m_how_many_x*m_how_many_y);
            std::copy(vals.begin(), vals.end(), m_values.begin());

            m_x_map_functor = XMapFunctor::deserialize(it_raw_data);
            m_y_map_functor = YMapFunctor::deserialize(it_raw_data);
            m_ix_map_functor = IXMapFunctor::deserialize(it_raw_data);
            m_iy_map_functor = IYMapFunctor::deserialize(it_raw_data);
        }

        /**
        * Operator ==
        *
        * @param[in] rhs a const reference to a 2D table of the same type
        *
        * @return true if rhs and *this are equal. false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator== (
            const generic_2d_table<
                RealType, VectorType,
                XMapFunctor, YMapFunctor,
                IXMapFunctor, IYMapFunctor> &rhs) const
        {
            return
                (m_how_many_x == rhs.m_how_many_x) &&
                (m_how_many_y == rhs.m_how_many_y) &&
                (m_values == rhs.m_values) &&
                (m_x_map_functor == rhs.m_x_map_functor) &&
                (m_y_map_functor == rhs.m_y_map_functor) &&
                (m_ix_map_functor == rhs.m_ix_map_functor) &&
                (m_iy_map_functor == rhs.m_iy_map_functor);
        }

        /**
        * Returns the i-th coodinate along x axis
        *
        * @param[in] i the index of the desired coordinate
        *
        * @return the i-th coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_coord(int i) const noexcept
        {
            return m_x_map_functor(i);
        }

        /**
        * Returns the j-th coodinate along y axis
        *
        * @param[in] j the index of the desired coordinate
        *
        * @return the j-th coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_coord(int j) const noexcept
        {
            return m_y_map_functor(j);
        }

        /**
        * Returns the function value at (i, j)
        *
        * @param[in] i the index of the desired value along x
        * @param[in] j the index of the desired value along y
        *
        * @return the value at (i,j)
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_val(int i, int j) const noexcept
        {
            return m_values[idx(i, j)];
        }

        /**
        * Returns the number of points along x
        *
        * @return the number of points along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int get_how_many_x() const noexcept
        {
            return m_how_many_x;
        }

        /**
        * Returns the minimum x coordinate
        *
        * @return the minimum x coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_min() const noexcept
        {
            return m_x_map_functor(0);
        }

        /**
        * Returns the maximum x coordinate
        *
        * @return the maximum x coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_max() const noexcept
        {
            return m_x_map_functor(m_how_many_x-1);
        }

        /**
        * Returns the size along x
        *
        * @return the size along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_x_size() const noexcept
        {
            return get_x_max() - get_x_min();
        }

        /**
        * Returns the size of the steps along x between i and i+1
        *
        * @return the size of the steps along x
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_dx(int i) const noexcept
        {
            return get_x_coord(i+1) - get_x_coord(i);
        }

        /**
        * Returns the number of points along y
        *
        * @return the number of points along y
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int get_how_many_y() const noexcept
        {
            return m_how_many_y;
        }

        /**
        * Returns the minimum y coordinate
        *
        * @return the minimum y coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_min() const noexcept
        {
            return m_y_map_functor(0);
        }

        /**
        * Returns the maximum y coordinate
        *
        * @return the maximum y coordinate
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_max() const noexcept
        {
            return m_y_map_functor(m_how_many_y-1);
        }

        /**
        * Returns the size along y
        *
        * @return the size along y
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_y_size() const noexcept
        {
            return get_y_max() - get_y_min();
        }

        /**
        * Returns the size of the steps along y
        *
        * @return the size of the steps along y
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType get_dy(int j) const noexcept
        {
            return get_y_coord(j+1) - get_y_coord(j);
        }

        /**
        * Returns an std::vector containing all the coordinates (not usable on
        * GPUs). Row major order is used.
        *
        * @return all the coordinates
        */
        std::vector<std::array<RealType,2>> get_all_coordinates() const noexcept
        {
            auto all_coords = std::vector<std::array<RealType,2>>(
                m_how_many_x*m_how_many_y, {0,0});
            int count = 0;
            for (int i = 0; i < m_how_many_x; ++i){
                for (int j = 0; j < m_how_many_y; ++j){
                    all_coords[count][0] = get_x_coord(i);
                    all_coords[count][1] = get_y_coord(j);
                    count++;
                }
            }
            return all_coords;
        }

        /**
        * Returns a const reference to the underlying Vector holding value data
        *
        * @return a const reference to the underlying Vector holding value data
        */
        const VectorType& get_values_reference() const noexcept
        {
            return m_values;
        }

        /**
        * Performs an interpolation a bilinear interpolation at given position.
        *
        * @param[in] where_x the position along x
        * @param[in] where_y the position along y
        * @return the result of the interpolation
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp(const RealType where_x, const RealType where_y) const noexcept
        {
            using namespace picsar::multi_physics::math;

            auto idx_x_left = m_ix_map_functor(where_x);
            if (idx_x_left == (m_how_many_x-1))
                idx_x_left = m_how_many_x-2;
            const auto idx_x_right = idx_x_left + 1;
            auto idx_y_left = m_iy_map_functor(where_y);
            if (idx_y_left == (m_how_many_y-1))
                idx_y_left = m_how_many_y-2;
            const auto idx_y_right = idx_y_left + 1;

            const auto xleft = get_x_coord(idx_x_left);
            const auto xright = get_x_coord(idx_x_right);
            const auto yleft = get_y_coord(idx_y_left);
            const auto yright = get_y_coord(idx_y_right);
            const auto f_xl_yl = m_values[idx(idx_x_left, idx_y_left)];
            const auto f_xl_yr = m_values[idx(idx_x_left, idx_y_right)];
            const auto f_xr_yl = m_values[idx(idx_x_right, idx_y_left)];
            const auto f_xr_yr = m_values[idx(idx_x_right, idx_y_right)];

            return utils::bilinear_interp(
                xleft, xright, yleft, yright,
                f_xl_yl, f_xl_yr, f_xr_yl, f_xr_yr,
                where_x, where_y);
        }

        /**
        * Performs an interpolation a linear interpolation at given position,
        * in the special case in which the second coordinate is exactly a grid
        * point along y
        *
        * @param[in] where_x the position along x
        * @param[in] j the index of the position along y
        * @return the result of the interpolation
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp_first_coord(RealType where_x, int j) const noexcept
        {
            using namespace picsar::multi_physics::math;

            auto idx_left = m_ix_map_functor(where_x);
            if (idx_left == (m_how_many_x-1))
                idx_left = m_how_many_x-2;
            const auto idx_right = idx_left + 1;

            const auto xleft = get_x_coord(idx_left);
            const auto xright = get_x_coord(idx_right);
            const auto left_val = m_values[idx(idx_left,j)];
            const auto right_val = m_values[idx(idx_right,j)];

            return utils::linear_interp(xleft, xright, left_val, right_val, where_x);
        }

        /**
        * Performs an interpolation a linear interpolation at given position,
        * in the special case in which the first coordinate is exactly a grid
        * point along x
        *
        * @param[in] i the position along x
        * @param[in] where_y the position along y
        * @return the result of the interpolation
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType interp_second_coord(int i, RealType where_y) const noexcept
        {
            using namespace picsar::multi_physics::math;

            auto idx_left = m_iy_map_functor(where_y);
            if (idx_left == (m_how_many_y-1))
                idx_left = m_how_many_y-2;
            const auto idx_right = idx_left + 1;

            const auto left_val = m_values[idx(i, idx_left)];
            const auto right_val = m_values[idx(i, idx_right)];
            const auto yleft = get_y_coord(idx_left);
            const auto yright = get_y_coord(idx_right);

            return utils::linear_interp(yleft, yright, left_val, right_val, where_y);
        }

        /**
        * Overwrites the (i,j) value (not designed for GPU usage)
        * Warning: no safety check on i and j is performed!
        *
        * @param[in] i index of the value to be overwritten
        * @param[in] j index of the value to be overwritten
        * @param[in] what value to be written at position (i,j)
        */
        PXRMP_FORCE_INLINE
        void set_val(const int i, const int j, const RealType what)
        {
            m_values[idx(i, j)] = what;
        }

        /**
        * Overwrites the i-th value (not designed for GPU usage)
        * Warning: no safety check on i is performed!
        *
        * @param[in] i index of the value to be overwritten
        * @param[in] what value to be written at position i
        */
        PXRMP_FORCE_INLINE
        void set_val(const int i, const RealType what)
        {
            m_values[i] = what;
        }

        /**
        * Returns a byte vector containing all the raw data of the
        * table (not usable on GPUs)
        *
        * @return a byte vector containing table data
        */
        std::vector<char> serialize() const
        {
            auto raw_data = std::vector<char>{};

            utils::serialization::put_in(
                static_cast<char>(sizeof(RealType)), raw_data);
            utils::serialization::put_in(m_how_many_x, raw_data);
            utils::serialization::put_in(m_how_many_y, raw_data);
            for (const auto val : m_values)
                utils::serialization::put_in(val, raw_data);
            for (const auto dd : m_x_map_functor.serialize())
                utils::serialization::put_in(dd, raw_data);
            for (const auto dd : m_y_map_functor.serialize())
                utils::serialization::put_in(dd, raw_data);
            for (const auto dd : m_ix_map_functor.serialize())
                utils::serialization::put_in(dd, raw_data);
            for (const auto dd : m_iy_map_functor.serialize())
                utils::serialization::put_in(dd, raw_data);

            return raw_data;
        }

    protected:

            int m_how_many_x = 0; /* how many grid points along x */
            int m_how_many_y = 0; /* how many grid points along y */
            VectorType m_values; /* values f(x,y) */
            XMapFunctor m_x_map_functor /* functor to map indices to coordinates on the first axis */;
            YMapFunctor m_y_map_functor /* functor to map indices to coordinates on the second axis */;
            IXMapFunctor m_ix_map_functor /* functor to map coordinates to indices on the first axis */;
            IYMapFunctor m_iy_map_functor/* functor to map coordinates to indices on the second axis */;

        /**
        * Function values are stored internally in a 1D vector.
        * This function calculates the index along this vector
        * from the coordinate indices i,j. Row major order is used.
        *
        * @param[in] i index along x
        * @param[in] j index along y
        * @return index along internal 1D vector
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int idx(int i, int j) const noexcept
        {
            return i*m_how_many_y + j;
        }
    };

    //__________________________________________________________________________


}
}
}



#endif //PICSAR_TABLES
