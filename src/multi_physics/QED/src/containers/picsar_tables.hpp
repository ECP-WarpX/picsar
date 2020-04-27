#ifndef PICSAR_EQUISPACED_TABLES
#define PICSAR_EQUISPACED_TABLES

#include<cmath>

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include "../utils/picsar_algo.hpp"

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
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        equispaced_1d_table(
            RealType x_min, RealType x_max, VectorType values):
            m_x_min{x_min}, m_x_max{x_max}, m_values{values}
            {
                m_how_many_x = static_cast<int>(values.size());
                m_x_size = x_max - x_min;
            };

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

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
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

    protected:
        RealType m_x_min;
        RealType m_x_max;
        RealType m_x_size;
        int m_how_many_x;
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

            const auto xleft = idx_left*m_x_size/(m_x_size-1) + m_x_min;
            const auto xright = idx_right*m_x_size/(m_x_size-1) + m_x_min;
            const auto left_val = m_values[idx(idx_left,j)];
            const auto right_val = m_values[idx(idx_right,j)];

            return linear_interp(xleft, xright, left_val, right_val, where_x);
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
            const auto yleft = idx_left*m_y_size/(m_y_size-1) + m_y_min;
            const auto yright = idx_right*m_y_size/(m_y_size-1) + m_y_min;

            return linear_interp(yleft, yright, left_val, right_val, where_y);
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


    protected:
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
