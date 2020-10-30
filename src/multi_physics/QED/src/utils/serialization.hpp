#ifndef PICSAR_MULTIPHYSICS_SERIALIZATION
#define PICSAR_MULTIPHYSICS_SERIALIZATION

//Should be included by all the src files of the library
#include "../qed_commons.h"

#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>
#include <cstring>

namespace picsar{
namespace multi_physics{
namespace utils{
namespace serialization{

    /**
    * This function transform a variable of type T into a vector of chars holding its
    * byte representation and it appends this vector at the end of an
    * existing vector of chars. T must be a POD type. This functions is not
    * usable on GPUs.
    *
    * @tparam T the variable type (must be POD)
    * @param[in, out] vec a reference to the vector to which the byte representation of val is appended
    */
    template<typename T>
    inline void put_in(T val, std::vector<char>& vec)
    {
        static_assert(std::is_pod<T>(), "Cannot serialize \
            non-POD types.");

        const auto* ptr_val = reinterpret_cast<char*>(&val);
        vec.insert(vec.end(), ptr_val, ptr_val+sizeof(T));
    }

    /**
    * This function extract a variable of type T from a byte vector, at the position
    * given by an iterator. The iterator is then advanced according to
    * the number of bytes read from the byte vector. T must be a POD type.
    * This functions is not usable on GPUs.
    *
    * @tparam T the variable type (must be POD)
    * @tparam CharIter the iterator type
    * @param[in, out] it the iterator to a byte vector
    * @return the variable extracted from the byte array
    */
    template<typename T, typename CharIter=std::vector<char>::const_iterator>
    inline T get_out(CharIter& it)
    {
        static_assert(std::is_pod<T>(), "Cannot extract \
            non-POD types from char vectors.");

        auto temp = std::array<char, sizeof(T)>{};
        std::copy(it, it + sizeof(T), temp.begin());
        it += sizeof(T);
        T res;
        std::memcpy(&res, temp.data(), sizeof(T));
        return res;
    }

    /**
    * This function extract several variables of type T from a byte vector,
    * at the position given by an iterator. The iterator is then advanced according to
    * the number of bytes read from the byte vector. T must be a POD type.
    * This functions is not usable on GPUs.
    *
    * @tparam T the variable type (must be POD)
    * @tparam CharIter the iterator type
    * @param[in, out] it the iterator to a byte vector
    * @param[in] how_many how many variables should be extracted
    * @return an std::vector<T> containing the extracted variables.
    */
    template<typename T, typename CharIter=std::vector<char>::const_iterator>
    inline std::vector<T> get_n_out(
        CharIter& it,
        int how_many)
    {
        static_assert(std::is_pod<T>(), "Cannot extract \
            non-POD types from char vectors.");

        auto res = std::vector<T>(how_many);
        std::generate(res.begin(), res.end(),
            [&](){
                auto temp = std::array<char, sizeof(T)>{};
                std::copy(it, it + sizeof(T), temp.begin());
                it += sizeof(T);
                T el;
                std::memcpy(&el, temp.data(), sizeof(T));
                return el;});
        return res;
    }

}
}
}
}

#endif //PICSAR_MULTIPHYSICS_SERIALIZATION
