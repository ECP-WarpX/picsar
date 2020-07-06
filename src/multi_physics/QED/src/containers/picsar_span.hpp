#ifndef PICSAR_MULTIPHYSICS_VIEW
#define PICSAR_MULTIPHYSICS_VIEW

//Should be included by all the src files of the library
#include "../qed_commons.h"

namespace picsar{
namespace multi_physics{
namespace containers{

    /**
    * This class implements a non-owning array
    *
    * @tparam T element type
    */
    template <typename T>
    class picsar_span
    {

    public:
        /**
        * Empty constructor.
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span(){}

        /**
        * Constructor requiring the size of the array and the pointer to raw data.
        *
        * @param[in] t_size size of the array
        * @param[in] ptr_data pointer to the first element of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_span(size_t t_size, T* ptr_data):
            m_size{t_size}, m_ptr_data{ptr_data}
        {}

       /**
        * Returns a reference to the i-th element of the array
        *
        * @param[in] i index of the desired element
        * @return a reference to the i-th element
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        T& operator [] (int i) noexcept
        {
            return m_ptr_data[i];
        }

        /**
        * Returs a const reference to the i-th element
        *
        * @param[in] i index of the desired element
        * @return a const reference to the i-th element
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const T& operator [] (int i) const noexcept
        {
            return m_ptr_data[i];
        }

       /**
        * Returs a const pointer to the underlying raw data array
        *
        * @return a const pointer to the underlying raw data array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* data() const noexcept
        {
            return m_ptr_data;
        }

       /**
        * Returs the size of the array
        *
        * @return the size of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr size_t size() const noexcept
        {
            return m_size;
        }

       /**
        * Returs a const pointer to the beginning of the array
        *
        * @return a const pointer to the first element of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* begin() const noexcept
        {
            return m_ptr_data;
        }

       /**
        * Returs a const pointer to the end of the array (i.e. the element after the last)
        *
        * @return a const pointer to the end of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr const T* end() const noexcept
        {
            return m_ptr_data+m_size;
        }

       /**
        * Returs a pointer to the beginning of the array
        *
        * @return a pointer to the first element of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* begin() noexcept
        {
            return m_ptr_data;
        }

       /**
        * Returs a pointer to the end of the array (i.e. the element after the last)
        *
        * @return a pointer to the end of the array
        */
        PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        constexpr T* end() noexcept
        {
            return m_ptr_data+m_size;
        }

        typedef T value_type;

        protected:
            size_t m_size = 0;  /* Array size */
            T* m_ptr_data = nullptr; /* Raw pointer to the array data */
        };
}
}
}

#endif //PICSAR_MULTIPHYSICS_VIEW
