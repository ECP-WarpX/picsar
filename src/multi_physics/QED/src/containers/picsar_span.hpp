#ifndef __PICSAR_MULTIPHYSICS_VIEW__
#define __PICSAR_MULTIPHYSICS_VIEW__

//This .hpp file contains the definition of a GPU-friendly STL-like array
//(thanks to Weiqun Zhang)

#include <cstddef>
#include <vector>

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
        picsar_span(){}


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span(size_t t_size, T* ptr_data):
            m_size{t_size}, m_ptr_data{ptr_data}
        {}


        template <size_t N>
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span(const picsar_array<T, N>& pic):
            m_size{N}, m_ptr_data{pic.data()}
        {}


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span(const std::vector<T>& vec):
            m_size{vec.size()}, m_ptr_data{vec.data()}
        {}


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T& operator [] (int i) const noexcept
        {
            return m_ptr_data[i];
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        T& operator [] (int i) noexcept
        {
            return m_ptr_data[i];
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T* data() const noexcept
        {
            return m_ptr_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        size_t size() const noexcept
        {
            return m_size;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T* begin() const noexcept
        {
            return m_ptr_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T* end() const noexcept
        {
            return m_ptr_data+m_size;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_FORCE_INLINE
        T* begin() noexcept
        {
            return m_ptr_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_FORCE_INLINE
        T* end() noexcept
        {
            return m_ptr_data;
        }


        private:
            T* m_ptr_data = nullptr;
            size_t m_size = 0;

        };
}
}
}

#endif //__PICSAR_MULTIPHYSICS_VIEW__
