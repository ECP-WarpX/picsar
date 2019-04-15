#ifndef __PICSAR_MULTIPHYSICS_VECTOR__
#define __PICSAR_MULTIPHYSICS_VECTOR__

//This .hpp file contains the definition of a GPU-friendly STL-like vector
//Dinamic resizing is NOT implemented (push_back, resize...).
//The (dangerous) possibility to construct using an existing raw pointer
//is provided

#include <cstddef>

//Should be included by all the src files of the library
#include "qed_commons.h"

//############################################### Declaration

namespace picsar{
    namespace multi_physics{

        template <typename T>
        class picsar_vector
        {

        public:
            PXRMP_GPU PXRMP_FORCE_INLINE
            picsar_vector(size_t size);

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

            PXRMP_GPU PXRMP_FORCE_INLINE
            const T& operator [] (int i) const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T& operator [] (int i) noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            T* data() const noexcept;

            PXRMP_GPU PXRMP_FORCE_INLINE
            size_t size() const noexcept;

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
picsar::multi_physics::picsar_vector&
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
picsar::multi_physics::picsar_vector&
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

#endif //__PICSAR_MULTIPHYSICS_ARRAY__
