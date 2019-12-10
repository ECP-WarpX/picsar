#ifndef __PICSAR_MULTIPHYSICS_ARRAY__
#define __PICSAR_MULTIPHYSICS_ARRAY__

//This .hpp file contains the definition of a GPU-friendly STL-like array
//(thanks to Weiqun Zhang)

#include <cstddef>

//Should be included by all the src files of the library
#include "../qed_commons.h"

#ifndef PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY
    #include <array>
#endif

//############################################### Declaration

namespace picsar{
namespace multi_physics{
namespace containers{

#ifdef PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY

    template <typename T, size_t N>
    class picsar_array
    {

    public:
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr picsar_array(){};


        // Constructor to allow initialization with a list {a,b,c,d,...}
        template<typename ... V>
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_array(V ... args):
        m_data{args...}
        {
            static_assert(sizeof...(args) == N, "Wrong number of elements");
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T& operator [] (int i) const noexcept
        {
            return m_data[i];
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        T& operator [] (int i) noexcept
        {
            return m_data[i];
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T* data() const noexcept
        {
            return m_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr size_t size() const noexcept
        {
            return N;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* begin() const noexcept
        {
            return m_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* end() const noexcept
        {
            return m_data+N;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* begin() noexcept
        {
            return m_data;
        }


        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* end() noexcept
        {
            return m_data+N;
        }


        private:
            T m_data[N];

        };
#else
    template <typename T, size_t N>
    using picsar_array = std::array<T, N>;
#endif

}
}
}

#endif //__PICSAR_MULTIPHYSICS_ARRAY__
