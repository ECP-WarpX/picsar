#ifndef PICSAR_MULTIPHYSICS_SERIALIZATION
#define PICSAR_MULTIPHYSICS_SERIALIZATION

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>

namespace picsar{
namespace multi_physics{
namespace utils{
namespace serialization{

    template<typename T>
    void put_in(T val, std::vector<char>& vec)
    {
        static_assert(std::is_pod<T>(), "Cannot serialize \
            non-POD types.");

        const auto* ptr_val = reinterpret_cast<char*>(&val);
        vec.insert(vec.end(), ptr_val, ptr_val+sizeof(T));
    }

    template<typename T, typename CharIter=std::vector<char>::const_iterator>
    T get_out(CharIter& it)
    {
        static_assert(std::is_pod<T>(), "Cannot extract \
            non-POD types from char vectors.");

        auto ptr_res = reinterpret_cast<const T*>(&(*it));
        it += sizeof(T);
        return *ptr_res;
    }

    template<typename T, typename CharIter=std::vector<char>::const_iterator>
    std::vector<T> get_n_out(
        CharIter& it,
        int how_many)
    {
        static_assert(std::is_pod<T>(), "Cannot extract \
            non-POD types from char vectors.");

        auto res = std::vector<T>(how_many);
        std::generate(res.begin(), res.end(),
            [&](){
                auto el = *reinterpret_cast<const T*>(&(*it));
                it += sizeof(T);
                return el;});
        return res;
    }



}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SERIALIZATION
