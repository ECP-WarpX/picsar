#ifndef PICSAR_MULTIPHYSICS_ARRAY
#define PICSAR_MULTIPHYSICS_ARRAY

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <cstddef>
//If GPUs are not used, picsar_arrays are just an alias for std::array
#ifndef PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY
    #include <array>
#endif

namespace picsar{
namespace multi_physics{
namespace containers{

#ifdef PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY

    /**
    * This class implements a GPU-friendly STL-like array.
    *
    * @tparam T type of the array name
    * @tparam N (size_t) the array size
    */
    template <typename T, size_t N>
    struct picsar_array
    {

        T m_data[N]; /* The underlying raw array */

       /**
        * Returs a const reference to the i-th element.
        *
        * @param[in] i the index of the element
        * @return a const reference to the i-th element
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T& operator [] (int i) const noexcept
        {
            return m_data[i];
        }

       /**
        * Returs a reference to the i-th element.
        *
        * @param[in] i the index of the element
        * @return a reference to the i-th element
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        T& operator [] (int i) noexcept
        {
            return m_data[i];
        }

       /**
        * Returs a const pointer to the underlying raw data array.
        *
        * @return a const pointer to the underlying raw data array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T* data() const noexcept
        {
            return m_data;
        }

       /**
        * Returs the size of the array.
        *
        * @return the size of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr int size() const noexcept
        {
            return static_cast<int>(N);
        }

       /**
        * Returs a const pointer to the beginning of the array.
        *
        * @return a const pointer to the first element of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* begin() const noexcept
        {
            return m_data;
        }

       /**
        * Returs a const pointer to the end of the array (i.e. the element after the last)
        *
        * @return a const pointer to the end of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* end() const noexcept
        {
            return m_data+N;
        }

       /**
        * Returs a pointer to the beginning of the array
        *
        * @return a pointer to the first element of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* begin() noexcept
        {
            return m_data;
        }

       /**
        * Returs a pointer to the end of the array (i.e. the element after the last)
        *
        * @return a pointer to the end of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* end() noexcept
        {
            return m_data+N;
        }

        };
#else

    // If PXRMP_INTERNAL_ENABLE_GPU_FRIENDLY_ARRAY is not defined,
    // picsar_array is just an alias for std::array
    template <typename T, size_t N>
    using picsar_array = std::array<T, N>;
#endif

}
}
}

#endif //PICSAR_MULTIPHYSICS_ARRAY
