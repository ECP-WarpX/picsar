#ifndef PICSAR_EQUISPACED_TABLES
#define PICSAR_EQUISPACED_TABLES

#include<cmath>

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include "../utils/picsar_algo.hpp"

namespace picsar{
namespace multi_physics{
namespace containers{

    template <typename RealType, class VectorType>
    class equispaced_1d_table{

    public:
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        equispaced_1d_table(
            RealType x_min, RealType x_max, VectorType values):
            m_x_min{x_min}, m_x_max{x_max}, m_values{values}
            {
                m_how_many_x = values.size();
                m_x_size = x_max - x_min;
            };

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType x_coord(const size_t i) const noexcept
        {
                return i*m_x_size/(m_how_many_x-1) + m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp(const RealType where_x) const noexcept
        {
            const auto idx_left = static_cast<size_t>(
                floor((m_how_many_x-1)*(where_x-m_x_min)/m_x_size));
            const auto idx_right = idx_left + 1;

            const auto xleft = x_coord(idx_left);
            const auto xright = x_coord(idx_right);
            const auto yleft = m_values[idx_left];
            const auto yright = m_values[idx_right];

            return utils::linear_interp(xleft, xright, yleft, yright, where_x);
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        void set_val(size_t where, RealType what)
        {
            m_values[where] = what;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        size_t get_how_many_x() const noexcept
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
        size_t m_how_many_x;
        VectorType m_values;
    };
/*
    template <typename RealType, class VectorType>
    class equispaced_2d_table
    {
    public:
        equispaced_2d_table(
            RealType x_min, RealType x_max,
            RealType y_min, RealType y_max,
            size_t how_many_x, size_t how_many_y,
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
        RealType x_coord(size_t i) const noexcept
        {
            return i*m_x_size/(m_how_many_x-1) + m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType y_coord(size_t j) const noexcept
        {
            return j*m_y_size/(m_how_many_y-1) + m_y_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp_first_coord(RealType where_x, size_t j) const
        {
            const auto idx_left = static_cast<size_t>(
                floor((m_x_size-1)*(where_x-m_x_min)/m_x_size));
            const auto idx_right = idx_left + 1;

            const auto xleft = idx_left*m_x_size/(m_x_size-1) + m_x_min;
            const auto xright = idx_right*m_x_size/(m_x_size-1) + m_x_min;
            const auto yleft = m_values[idx(idx_left,j)];
            const auto yright = m_values[idx(idx_right,j)];

            return linear_interp(xleft, xright, yleft, yright, where_x);
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType interp_second_coord(size_t i, RealType where_y) const
        {
            const auto idx_left = static_cast<size_t>(
                floor((m_y_size-1)*(where_y-m_y_min)/m_y_size));
            const auto idx_right = idx_left + 1;

            const auto xleft = m_values[idx(i, idx_left)];
            const auto xright = m_values[idx(i, idx_right)];
            const auto yleft = idx_left*m_y_size/(m_y_size-1) + m_y_min;
            const auto yright = idx_right*m_y_size/(m_y_size-1) + m_y_min;

            return linear_interp(xleft, xright, yleft, yright, where_x);
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        void set_val(size_t where_x, size_t where_y, RealType what)
        {
            m_values[idx(where_x, where_y)] = what;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        size_t get_how_many_x(size_t i) const
        {
            return m_how_many_x;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_min() const
        {
            return m_x_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_x_max() const
        {
            return m_x_max;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        size_t get_how_many_y(size_t i) const
        {
            return m_how_many_y;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_y_min() const
        {
            return m_y_min;
        }

        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        RealType get_y_max() const
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
            size_t m_how_many_x;
            size_t m_how_many_y;
            VectorType m_values;

            PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            size_t idx(size_t i, size_t j) const
            {
                return i*m_y_size + j;
            }
    };

    */
}
}
}



#endif //PICSAR_TABLES
