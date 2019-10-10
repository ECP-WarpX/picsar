#ifndef __PICSAR_MULTIPHYSICS_ARRAY__
#define __PICSAR_MULTIPHYSICS_ARRAY__

//This .hpp file contains the definition of a GPU-friendly STL-like array
//(thanks to Weiqun Zhang)

#include <cstddef>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        template <typename T, size_t N>
        class picsar_array
        {

        public:
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_array();

            // Constructor to allow initialization with a list {a,b,c,d,...}
            template<typename ... V>
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_array(V ... args);

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T& operator [] (int i) const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T& operator [] (int i) noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T* data() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            size_t size() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T* begin() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T* end() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T* begin() noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T* end() noexcept;


        private:
            T arr[N];


        };

    }
}

//############################################### Implementation

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_array<T, N>::picsar_array()
{}

// Constructor to allow initialization with a list {a,b,c,d,...}
template <typename T, size_t N>
template<typename ... V>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_array<T, N>::picsar_array(V ... args):
arr{args...}
{
    static_assert(sizeof...(args) == N, "Wrong number of elements");
}

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
const T&
picsar::multi_physics::picsar_array<T, N>::operator [](int i) const noexcept
{
    return arr[i];
}


template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
T&
picsar::multi_physics::picsar_array<T, N>::operator [](int i) noexcept
{
    return arr[i];
}


template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
const T*
picsar::multi_physics::picsar_array<T, N>:: data() const noexcept
{
    return arr;
}

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
size_t
picsar::multi_physics::picsar_array<T, N>::size() const noexcept
{
    return N;
}

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
const T*
picsar::multi_physics::picsar_array<T, N>:: begin() const noexcept
{
    return arr;
}

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
T*
picsar::multi_physics::picsar_array<T, N>:: begin() noexcept
{
    return arr;
}

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
const T*
picsar::multi_physics::picsar_array<T, N>:: end() const noexcept
{
    return arr+N;
}

template <typename T, size_t N>
PXRMP_GPU PXRMP_FORCE_INLINE
T*
picsar::multi_physics::picsar_array<T, N>:: end() noexcept
{
    return arr+N;
}

#endif //__PICSAR_MULTIPHYSICS_ARRAY__
