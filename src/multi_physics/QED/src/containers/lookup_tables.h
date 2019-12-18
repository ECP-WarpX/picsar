#ifndef PICSAR_MULTIPHYSICS_LOOKUP_TABLE
#define PICSAR_MULTIPHYSICS_LOOKUP_TABLE

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include "picsar_array.hpp"

namespace picsar{
namespace multi_physics{
namespace containers{

    template <class ContainerType, size_t N>
    struct lookup_table{
        picsar_array<ContainerType, N> coords;
        picsar_array<size_t, N> data_sizes;
        ContainerType data;
    };

}
}
}

#endif //PICSAR_MULTIPHYSICS_LOOKUP_TABLE