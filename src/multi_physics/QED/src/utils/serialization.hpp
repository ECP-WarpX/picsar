#ifndef PICSAR_MULTIPHYSICS_SERIALIZATION
#define PICSAR_MULTIPHYSICS_SERIALIZATION

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <vector>
#include <array>
#include <algorithm>

namespace picsar{
namespace multi_physics{
namespace utils{
namespace serialization{

    template<typename T>
    void put_in(T val, std::vector<char>& vec)
    {
        const auto* ptr_val = reinterpret_cast<char*>(&val);
        vec.insert(vec.end(), ptr_val, ptr_val+sizeof(T));
    }

    template<typename T, typename CharIter=std::vector<char>::const_iterator>
    T get_out(CharIter& it)
    {
        auto ptr_res = reinterpret_cast<const T*>(&(*it));
        it += sizeof(T);
        return *ptr_res;
    }

    template<typename T, typename CharIter=std::vector<char>::const_iterator>
    std::vector<T> get_n_out(
        CharIter& it,
        int how_many)
    {
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
