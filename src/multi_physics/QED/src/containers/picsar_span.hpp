#ifndef __PICSAR_MULTIPHYSICS_VIEW__
#define __PICSAR_MULTIPHYSICS_VIEW__

#include <cstddef>
#include <type_traits>

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include "picsar_array.hpp"

//############################################### Declaration

namespace picsar{
namespace multi_physics{
namespace containers{

    template <typename T>
    class picsar_span
    {

    public:
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span():
            m_size{0}, m_ptr_data{nullptr}{}


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span(size_t t_size, T* ptr_data):
            m_size{t_size}, m_ptr_data{ptr_data}
        {}


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        T& operator [] (int i) const noexcept
        {
            return m_ptr_data[i];
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* data() const noexcept
        {
            return m_ptr_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr size_t size() const noexcept
        {
            return m_size;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* begin() const noexcept
        {
            return m_ptr_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* end() const noexcept
        {
            return m_ptr_data+m_size;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* begin() noexcept
        {
            return m_ptr_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* end() noexcept
        {
            return m_ptr_data;
        }


        protected:
            T* m_ptr_data;
            size_t m_size;

        };
}
}
}

#endif //__PICSAR_MULTIPHYSICS_VIEW__
