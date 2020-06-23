#ifndef PICSAR_EQUISPACED_TABLES
#define PICSAR_EQUISPACED_TABLES

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include "../utils/picsar_algo.hpp"

#include "../utils/serialization.hpp"

#include "../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace containers{

    /**
    * This class implements a generic equispace 1D lookup table
    *
    * @tparam RealType the floating point type to be used (e.g. double or float)
    * @tparam VectorType the vector type to be used (e.g. std::vector<double>)
    */
    template <typename RealType, class VectorType>
    class equispaced_1d_table{

    public:

        /**
        * Constructor
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
            };

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        equispaced_1d_table(){};

        equispaced_1d_table(const std::vector<char>& raw_data)
        {
            using namespace picsar::multi_physics::utils;

            constexpr size_t min_size =
                sizeof(char)+//magic_number
                sizeof(m_x_min)+sizeof(m_x_max)+//xmin, xmax
                sizeof(m_how_many_x);// m_how_many_x

            if (raw_data.size() < min_size)
                throw "Binary data is too small to be a 1D table.";

            auto it_raw_data = raw_data.begin();

            if (serialization::get_out<char>(it_raw_data) !=
                static_cast<char>(sizeof(RealType))){
                throw "Mismatch between RealType used to write and to read the 1D table";
            }

            m_x_min = serialization::get_out<RealType>(it_raw_data);
            m_x_max = serialization::get_out<RealType>(it_raw_data);
            m_x_size = m_x_max - m_x_min;
            if(m_x_size < 0)
                throw "raw_data contains invalid data.";

            m_how_many_x = serialization::get_out<int>(it_raw_data);
            if(m_how_many_x <= 0)
                throw "raw_data contains invalid data.";
            m_values = VectorType(m_how_many_x);
            auto vals =
                serialization::get_n_out<RealType>(it_raw_data, m_how_many_x);
            std::copy(vals.begin(), vals.end(), m_values.begin());
        };

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        bool operator== (
            const equispaced_1d_table<RealType, VectorType> &b) const
        {
            return
                (m_x_min == b.m_x_min) &&
                (m_x_max == b.m_x_max) &&
                (m_x_size == b.m_x_size) &&
                (m_how_many_x == b.m_how_many_x) &&
                (m_values == b.m_values);
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_coord(const int i) const noexcept
        {
                return i*m_x_size/(m_how_many_x-1) + m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_val(const size_t i) const noexcept
        {
                return m_values[i];
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp(const RealType where_x) const noexcept
        {
            const auto idx_left = static_cast<int>(
                floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            if (idx_left == (m_how_many_x-1))
                return  m_values[m_how_many_x-1];
            const auto idx_right = idx_left + 1;

            const auto xleft = get_x_coord(idx_left);
            const auto xright = get_x_coord(idx_right);
            const auto yleft = m_values[idx_left];
            const auto yright = m_values[idx_right];

            return utils::linear_interp(xleft, xright, yleft, yright, where_x);
        }

        void set_val(int i, RealType what)
        {
            m_values[i] = what;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        int get_how_many_x() const noexcept
        {
            return m_how_many_x;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_min() const noexcept
        {
            return m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_max() const noexcept
        {
            return m_x_max;
        }

        std::vector<RealType> get_all_coordinates()  const noexcept
        {
            auto all_coords = std::vector<RealType>(m_how_many_x);
            for (int i = 0; i < m_how_many_x; ++i){
                all_coords[i] = get_x_coord(i);
            }
            return all_coords;
        }

        std::vector<char> serialize() const
        {
            auto raw_data = std::vector<char>{};

            utils::serialization::put_in(
                static_cast<char>(sizeof(RealType)), raw_data);
            utils::serialization::put_in(m_x_min, raw_data);
            utils::serialization::put_in(m_x_max, raw_data);
            utils::serialization::put_in(m_how_many_x, raw_data);
            for (auto val : m_values)
                utils::serialization::put_in(val, raw_data);

            return raw_data;
        }

        RealType m_x_min = 0.0;
        RealType m_x_max = 0.0;
        RealType m_x_size = 0.0;
        int m_how_many_x = 0;
        VectorType m_values;
    };

    template <typename RealType, class VectorType>
    class equispaced_2d_table
    {
    public:
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
            };

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_coord(int i) const noexcept
        {
            return i*m_x_size/(m_how_many_x-1) + m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_y_coord(int j) const noexcept
        {
            return j*m_y_size/(m_how_many_y-1) + m_y_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_val(int i, int j) const noexcept
        {
            return m_values[idx(i, j)];
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp_first_coord(RealType where_x, int j) const noexcept
        {
            auto idx_left = static_cast<int>(
                floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            if (idx_left == (m_how_many_x-1))
                idx_left = m_how_many_x-2;
            const auto idx_right = idx_left + 1;

            const auto xleft = (idx_left*m_x_size)/(m_how_many_x-1) + m_x_min;
            const auto xright = (idx_right*m_x_size)/(m_how_many_x-1) + m_x_min;
            const auto left_val = m_values[idx(idx_left,j)];
            const auto right_val = m_values[idx(idx_right,j)];

            return utils::linear_interp(xleft, xright, left_val, right_val, where_x);
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp_second_coord(int i, RealType where_y) const noexcept
        {
            auto idx_left = static_cast<int>(
                floor((m_how_many_y-1)*(where_y-m_y_min)/m_y_size));
            if (idx_left == (m_how_many_y-1))
                idx_left = m_how_many_y-2;
            const auto idx_right = idx_left + 1;

            const auto left_val = m_values[idx(i, idx_left)];
            const auto right_val = m_values[idx(i, idx_right)];
            const auto yleft = (idx_left*m_y_size)/(m_how_many_y-1) + m_y_min;
            const auto yright = (idx_right*m_y_size)/(m_how_many_y-1) + m_y_min;

            return utils::linear_interp(yleft, yright, left_val, right_val, where_y);
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp(const RealType where_x, const RealType where_y) const noexcept
        {
            auto idx_x_left = static_cast<int>(
                floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            if (idx_x_left == (m_how_many_x-1))
                idx_x_left = m_how_many_x-2;
            const auto idx_x_right = idx_x_left + 1;
            auto idx_y_left = static_cast<int>(
                floor((m_how_many_y-1)*(where_y-m_y_min)/m_y_size));
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

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        void set_val(const int i, const int j, const RealType what)
        {
            m_values[idx(i, j)] = what;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        int get_how_many_x() const noexcept
        {
            return m_how_many_x;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_min() const noexcept
        {
            return m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_max() const noexcept
        {
            return m_x_max;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        int get_how_many_y() const noexcept
        {
            return m_how_many_y;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_y_min() const noexcept
        {
            return m_y_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_y_max() const noexcept
        {
            return m_y_max;
        }

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

            RealType m_x_min;
            RealType m_x_max;
            RealType m_y_min;
            RealType m_y_max;
            RealType m_x_size;
            RealType m_y_size;
            int m_how_many_x;
            int m_how_many_y;
            VectorType m_values;

            PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            int idx(int i, int j) const noexcept
            {
                return i*m_how_many_y + j;
            }
    };
}
}
}



#endif //PICSAR_TABLES
