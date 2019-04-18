#ifndef __PICSAR_MULTIPHYSICS_VECTOR__
#define __PICSAR_MULTIPHYSICS_VECTOR__

//This .hpp file contains the definition of a GPU-friendly STL-like vector
//Resizing is implemented only for CPU.
//The (dangerous) possibility to construct using an existing raw pointer
//is provided

#include <cstddef>
#include <vector>
#include <algorithm>
#include <stdexcept>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        template <typename T>
        class picsar_vector
        {
        public:
            typedef T value_type;

            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector(size_t size);

            // Constructor to allow initialization with an STL vector (not for GPU)
            picsar_vector(const std::vector<T>& vec);

            //Empty initialization (not for GPU)
            picsar_vector();

            // Constructor to allow initialization with a list {a,b,c,d,...}
            //(not for GPU)
            template<typename ... V>
            picsar_vector(V ... args);

            // Constructor to allow initialization with a raw pointer
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector(size_t size, T* raw);

            //Copy constructor  (a deep copy)
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector(const picsar_vector& other);

            //Move constructor
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector(picsar_vector&& other);

            //Assignment operator  (a deep copy)
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector&  operator= (const picsar_vector& other);

            //Move assignment operator
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector&  operator= (picsar_vector&& other);

            PXRMP_GPU PXRMP_FORCE_INLINE
            ~picsar_vector();

            //Resize operation (not for GPU)
            void resize(size_t new_size);

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T& operator [] (int i) const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T& operator [] (int i) noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T* data() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            size_t size() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T* begin();

            PXRMP_GPU PXRMP_FORCE_INLINE
            T* end();

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T* begin() const;

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T* end() const;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T& back();

            PXRMP_GPU PXRMP_FORCE_INLINE
            T& front();

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T& back() const;

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T& front() const;


        private:
            T* v_data;
            size_t v_size;
            //If the class has been instantiated using a raw pointer, it should
            //*NOT* delete v_data in the destructor
            bool should_manage_memory;

        };

    }
}

//############################################### Implementation

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>::picsar_vector(size_t size)
{
    v_size = size;
    v_data = new T[size];

    should_manage_memory = true;
}


// Constructor to allow initialization with an STL vector (not for GPU)
template <typename T>
picsar::multi_physics::picsar_vector<T>::picsar_vector(const std::vector<T>& vec)
{
    v_size = vec.size();
    v_data = new T[v_size];

    std::copy(vec.begin(), vec.end(), v_data);

    should_manage_memory = true;
}

//Empty initialization (not for GPU)
template <typename T>
picsar::multi_physics::picsar_vector<T>::picsar_vector()
{
    v_size = 0;
    v_data = nullptr;

    should_manage_memory = true;
}

// Constructor to allow initialization with a list {a,b,c,d,...}
template <typename T>
template<typename ... V>
picsar::multi_physics::picsar_vector<T>::picsar_vector(V ... args):
v_size{sizeof...(args)}
{
    v_data = new T[v_size];
    std::vector<T> ttemp{args...};

    std::copy(ttemp.begin(), ttemp.end(), v_data);

    should_manage_memory = true;
}

// Constructor to allow initialization with a raw pointer
// CAUTION!!! This function is dangerous, since the user is finally responsible
// to allocate and free memory. It is advisable to use this method only inside
// CUDA kernels.
template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>::picsar_vector(size_t size, T* raw)
{
    v_size = size;
    v_data = raw;

    should_manage_memory = false;
}


//Copy constructor (a deep copy)
template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>::
picsar_vector(const picsar_vector& other)
{
    v_size = other.v_size;

    v_data = new T[other.v_size];
    for(size_t i = 0; i < v_size; i++)
        v_data[i] = other.v_data[i];

    should_manage_memory = true;
}


//Move constructor
template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>::picsar_vector(picsar_vector&& other)
{
    v_size = other.v_size;
    v_data = other.v_data;
    should_manage_memory = other.should_manage_memory;

    other.v_size = 0;
    other.v_data = nullptr;
    other.should_manage_memory = false;
}


//Assignment operator (a deep copy)
template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>&
picsar::multi_physics::picsar_vector<T>::operator= (const picsar_vector& other)
{
    v_size = other.v_size;

    v_data = new T[other.v_size];
    for(size_t i = 0; i < v_size; i++)
        v_data[i] = other.v_data[i];

    should_manage_memory = true;

    return *this;
}


//Move assignment operator
template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>&
picsar::multi_physics::picsar_vector<T>::operator= (picsar_vector&& other)
{
    v_size = other.v_size;
    v_data = other.v_data;
    should_manage_memory = other.should_manage_memory;

    other.v_size = 0;
    other.v_data = nullptr;
    other.should_manage_memory = false;

    return *this;
}


//If the class has been instantiated using a raw pointer, it should
//*NOT* delete v_data in the destructor
template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
picsar::multi_physics::picsar_vector<T>::~picsar_vector()
{
    if(should_manage_memory)
        delete[] v_data;
}

//Resize operation (not for GPU)
template <typename T>
void picsar::multi_physics::picsar_vector<T>::resize(size_t new_size)
{
    if(!should_manage_memory)
        throw std::logic_error("Cannot resize an underlying non owning pointer");

    T* temp_data = new T[new_size];

    size_t how_many_copy = size;

    if(size > new_size)
        how_many_copy = new_size;

    std::copy(v_data, v_data+how_many_copy, temp_data);
    std::swap(v_data, temp_data);

    delete[] temp_data;

    size = new_size;
}


template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
const T&
picsar::multi_physics::picsar_vector<T>::operator[] (int i) const noexcept
{
    return v_data[i];
}


template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
T& picsar::multi_physics::picsar_vector<T>::operator[] (int i) noexcept
{
    return v_data[i];
}


template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
T* picsar::multi_physics::picsar_vector<T>::data() const noexcept
{
    return v_data;
}


template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
size_t
picsar::multi_physics::picsar_vector<T>::size() const noexcept
{
    return v_size;
}

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
T*
picsar::multi_physics::picsar_vector<T>::begin()
{
    return v_data;
}

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
T*
picsar::multi_physics::picsar_vector<T>::end()
{
    return v_data + v_size;
}


template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
const T*
picsar::multi_physics::picsar_vector<T>::begin() const
{
    return v_data;
}

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
const T*
picsar::multi_physics::picsar_vector<T>::end() const
{
    return v_data + v_size;
}

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
T&
picsar::multi_physics::picsar_vector<T>::back()
{
    return v_data[0];
}

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
T&
picsar::multi_physics::picsar_vector<T>::front()
{
    return v_data[v_size-1];
}


template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
const T&
picsar::multi_physics::picsar_vector<T>::back() const
{
    return v_data[0];
}

template <typename T>
PXRMP_GPU PXRMP_FORCE_INLINE
const T&
picsar::multi_physics::picsar_vector<T>::front() const
{
    return v_data[v_size-1];
}

#endif //__PICSAR_MULTIPHYSICS_ARRAY__
