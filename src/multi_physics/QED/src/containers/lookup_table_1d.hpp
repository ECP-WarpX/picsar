#ifndef PICSAR_MULTIPHYSICS_LOOKUP_TABLE_1D
#define PICSAR_MULTIPHYSICS_LOOKUP_TABLE_1D

//Should be included by all the src files of the library
#include "../qed_commons.h"

namespace picsar{
namespace multi_physics{
namespace containers{

    template<typename T>
    struct picsar_pair
    {
        T first;
        T second;
    };

    template <class ContainerType>
    class lookup_table_1d{
    public:
        typedef typename ContainerType::value_type value_type;

        PXRMP_INTERNAL_GPU_DECORATOR
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        lookup_table_1d(ContainerType coords, ContainerType data):
            m_coords{coords}, m_data{data}{};

        PXRMP_INTERNAL_GPU_DECORATOR
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_pair<value_type> get_coord_val_at(int i) const noexcept
        {
            return picsar_pair<value_type>{m_coords[i], m_data[i]};
        }

        PXRMP_INTERNAL_GPU_DECORATOR
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        picsar_pair<value_type> set_coord_val_at(
            picsar_pair<value_type>& pair, int i) noexcept
        {
            m_coords[i] = pair.first;
            m_data[i] = pair.second;
        }

        PXRMP_INTERNAL_GPU_DECORATOR
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const value_type* get_coord_ptr() const noexcept
        {
            return m_coords.data_ptr;
        }

        PXRMP_INTERNAL_GPU_DECORATOR
        PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
        const value_type* get_data_ptr() const noexcept
        {
            return m_data.data_ptr;
        }

    private:
        ContainerType m_coords;
        ContainerType m_data;
    };

}
}
}

#endif