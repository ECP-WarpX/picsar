#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_DETAIL
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_DETAIL

//This .hpp file contains auxiliary functors used for the "tail-optimized"
//Quantum Synchrotron lookup table

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Uses serialization
#include "picsar_qed/utils/serialization.hpp"
//Uses log and exp
#include "picsar_qed/math/cmath_overloads.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{
namespace detail{

    /**
    * This class implements a linear functor to be used in tail-optimized lookup tables.
    * The functor maps integers between 0 and zsize-1 to linearly spaced values between
    * zmin and zmax. It also contains methods to serialize and deserialize itself.
    *
    * @tparam RealType the floating point type to be used
    */
    template <typename RealType>
    class LinFunctor
    {
    public:

        /**
        * Empty constructor
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        LinFunctor(){}

        /**
        * Constructor
        *
        * @param[in] zsize the number of points
        * @param[in] zmin the lower extremum
        * @param[in] zmax the upper extremum
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        LinFunctor(int zsize, RealType zmin, RealType zmax) :
            m_zmin{zmin}
        {
            m_coeff = (zmax - zmin)/(zsize - 1);
        }

        /**
        * Operator()
        *
        * @param[in] i an integer value
        *
        * @return zmin + i*(zmax-zmin)/(zsize-1)
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType operator()(const int i) const
        {
            return m_zmin + i * m_coeff;
        }

        /**
        * Operator==
        *
        * @param[in] rhs a const reference to a LinFunctor<RealType>
        *
        * @return true if rhs and the current functor are equal, false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator==(const LinFunctor<RealType> &rhs) const
        {
            return (this->m_zmin == rhs.m_zmin) &&
                   (this->m_coeff == rhs.m_coeff);
        }

        /**
        * Serializes the functor (not suitable for GPUs)
        *
        * @return an std::vector<char> containing the binary representation of the functor
        */
        std::vector<char> serialize() const
        {
            using namespace utils;

            std::vector<char> raw;

            serialization::put_in(m_zmin, raw);
            serialization::put_in(m_coeff, raw);

            return raw;
        }

        /**
        * Deserializes a LinFunctor object from binary data (not suitable for GPUs)
        *
        * @tparam CharIter a char interator type
        *
        * @param[in,out] it a char interator
        *
        * @return a LinFunctor<RealType> object
        */
        template <typename CharIter>
        static LinFunctor<RealType> deserialize(CharIter &it)
        {
            using namespace utils;

            auto functor = LinFunctor<RealType>{};

            functor.m_zmin = serialization::get_out<RealType>(it);
            functor.m_coeff = serialization::get_out<RealType>(it);

            return functor;
        }

    private:

        RealType m_zmin;
        RealType m_coeff;
    };


    /**
    * This class implements an inverse linear functor to be used in tail-optimized lookup tables.
    * The functor is essentially the inverse of LinFunctor. It also contains methods
    * to serialize and deserialize itself.
    *
    * @tparam RealType the floating point type to be used
    */
    template <typename RealType>
    class ILinFunctor
    {
    public:

        /**
        * Empty constructor
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        ILinFunctor(){}

        /**
        * Constructor
        *
        * @param[in] zsize the number of points
        * @param[in] zmin the lower extremum
        * @param[in] zmax the upper extremum
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        ILinFunctor(int zsize, RealType zmin, RealType zmax) :
            m_zmin{zmin}
        {
            m_coeff = (zsize - 1)/(zmax - zmin);
        }

        /**
        * Operator()
        *
        * @param[in] z a real value
        *
        * @return floor(z - zmin)*(zsize - 1)/(zmax - zmin)
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType operator()(const RealType z) const
        {
            using namespace picsar::multi_physics::math;

            return static_cast<int>(
                math::m_floor((z - m_zmin)*m_coeff));
        }

        /**
        * Operator==
        *
        * @param[in] rhs a const reference to an ILinFunctor<RealType>
        *
        * @return true if rhs and the current functor are equal, false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator==(const ILinFunctor<RealType> &rhs) const
        {
            return (this->m_zmin == rhs.m_zmin) &&
                   (this->m_coeff == rhs.m_coeff);
        }

        /**
        * Serializes the functor (not suitable for GPUs)
        *
        * @return an std::vector<char> containing the binary representation of the functor
        */
        std::vector<char> serialize() const
        {
            using namespace utils;

            std::vector<char> raw;

            serialization::put_in(m_zmin, raw);
            serialization::put_in(m_coeff, raw);

            return raw;
        }

        /**
        * Deserializes a ILinFunctor object from binary data (not suitable for GPUs)
        *
        * @tparam CharIter a char interator
        *
        * @param[in,out] it a char interator
        *
        * @return a ILinFunctor<RealType> object
        */
        template <typename CharIter>
        static ILinFunctor<RealType> deserialize(CharIter &it)
        {
            using namespace utils;

            auto functor = ILinFunctor<RealType>{};

            functor.m_zmin = serialization::get_out<RealType>(it);
            functor.m_coeff = serialization::get_out<RealType>(it);

            return functor;
        }

    private:

        RealType m_zmin;
        RealType m_coeff;
    };


    /**
    * This class implements a piecewise functor to be used in tail-optimized lookup tables.
    * The functor maps integers between 0 and zfirst-1 to linearly spaced values between
    * zmin and zswitch. It maps integers between zswitch and zsize-1 to logarithmically
    * spaced points up to zmax. It also contains methods to serialize and deserialize itself.
    *
    * @tparam RealType the floating point type to be used
    */
    template <typename RealType>
    class TailOptFunctor
    {
    public:

        /**
        * Empty constructor
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        TailOptFunctor(){}

        /**
        * Constructor
        *
        * @param[in] zsize the total number of points
        * @param[in] zfirst the number of points to be used in the first interval
        * @param[in] zmin the lower extremum
        * @param[in] zmax the upper extremum
        * @param[in] zswitch the z value corresponding to the switch between the two pieces of the function
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        TailOptFunctor(const int zsize, const int zfirst,
            const RealType zmin, const RealType zmax, const RealType zswitch) :
                m_zsize{zsize}, m_zfirst{zfirst}, m_zmin{zmin}
        {
            m_exp_zswitch = math::m_exp(zswitch);
            m_coeff_first = (zswitch - zmin) / (zfirst - 1);
            m_coeff_second = (math::m_exp(zmax) - m_exp_zswitch) / (zsize - zfirst);
        }

        /**
        * Operator()
        *
        * @param[in] i an integer value
        *
        * @return if i < zfirst,  zmin + i*(zswitch-zmin)/(zfirst-1), otherwise log( exp(zswitch) + (i+1-zfirst)*(exp(zmax) - exp(zswitch))/(zsize-zfirst) )
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType operator()(const int i) const
        {
            using namespace picsar::multi_physics::math;

            if (i < m_zfirst)
            {
                return m_zmin + i * m_coeff_first;
            }
            else
            {
                return math::m_log(
                    m_exp_zswitch + ((i + 1) - m_zfirst) * m_coeff_second);
            }
        }

        /**
        * Operator==
        *
        * @param[in] rhs a const reference to a TailOptFunctor<RealType>
        *
        * @return true if rhs and the current functor are equal, false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator==(const TailOptFunctor<RealType> &rhs) const
        {
            return (this->m_zsize == rhs.m_zsize) &&
                   (this->m_zfirst == rhs.m_zfirst) &&
                   (this->m_zmin == rhs.m_zmin) &&
                   (this->m_exp_zswitch == rhs.m_exp_zswitch) &&
                   (this->m_coeff_first == rhs.m_coeff_first) &&
                   (this->m_coeff_second == rhs.m_coeff_second);
        }

        /**
        * Serializes the functor (not suitable for GPUs)
        *
        * @return an std::vector<char> containing the binary representation of the functor
        */
        std::vector<char> serialize() const
        {
            using namespace utils;

            std::vector<char> raw;

            serialization::put_in(m_zsize, raw);
            serialization::put_in(m_zfirst, raw);
            serialization::put_in(m_zmin, raw);
            serialization::put_in(m_exp_zswitch, raw);
            serialization::put_in(m_coeff_first, raw);
            serialization::put_in(m_coeff_second, raw);

            return raw;
        }

        /**
        * Deserializes a TailOptFunctor object from binary data (not suitable for GPUs)
        *
        * @tparam CharIter a char interator type
        *
        * @param[in,out] it a char interator
        *
        * @return a LinFunctor<RealType> object
        */
        template <typename CharIter>
        static TailOptFunctor<RealType> deserialize(CharIter &it)
        {
            using namespace utils;

            auto functor = TailOptFunctor<RealType>{};

            functor.m_zsize = serialization::get_out<int>(it);
            functor.m_zfirst = serialization::get_out<int>(it);
            functor.m_zmin = serialization::get_out<RealType>(it);
            functor.m_exp_zswitch = serialization::get_out<RealType>(it);
            functor.m_coeff_first = serialization::get_out<RealType>(it);
            functor.m_coeff_second = serialization::get_out<RealType>(it);

            return functor;
        }

    private:

        int m_zsize;
        int m_zfirst;
        RealType m_zmin;
        RealType m_exp_zswitch;
        RealType m_coeff_first;
        RealType m_coeff_second;
    };

    /**
    * This class essentially implements the inverse functor of TailOptFunctor.
    * It also contains methods to serialize and deserialize itself.
    *
    * @tparam RealType the floating point type to be used
    */
    template <typename RealType>
    class ITailOptFunctor
    {
    public:

        /**
        * Empty constructor
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        ITailOptFunctor(){}

        /**
        * Constructor
        *
        * @param[in] zsize the total number of points
        * @param[in] zfirst the number of points to be used in the first interval
        * @param[in] zmin the lower extremum
        * @param[in] zmax the upper extremum
        * @param[in] zswitch the z value corresponding to the switch between the two pieces of the function
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        ITailOptFunctor(const int zsize, const int zfirst,
            const RealType zmin, const RealType zmax, const RealType zswitch) :
                m_zsize{zsize}, m_zfirst{zfirst}, m_zmin{zmin}, m_zswitch{zswitch}
        {
            m_exp_zswitch = math::m_exp(zswitch);
            m_coeff_first = (zfirst - 1) / (zswitch - zmin);
            m_coeff_second = (zsize - zfirst) / (math::m_exp(zmax) - m_exp_zswitch);
        }

        /**
        * Operator()
        *
        * @param[in] i an integer value
        *
        * @return if z < zswitch,  floor((z-zmin)*(zfirst-1))/(zswitch-zmin)), otherwise floor((exp(z)-exp(zswitch))*(zsize-zfirst)/(exp(zmax)-exp(zswitch)) + (zfirst-1))
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int operator()(const RealType z) const
        {
            if (z < m_zswitch)
            {
                return static_cast<int>(math::m_floor(
                    (z - m_zmin) * m_coeff_first));
            }
            else
            {
                const auto exp_z = math::m_exp(z);
                return static_cast<int>(math::m_floor(
                    (exp_z - m_exp_zswitch) * m_coeff_second + (m_zfirst - 1)));
            }
        }

        /**
        * Operator==
        *
        * @param[in] rhs a const reference to an ITailOptFunctor<RealType>
        *
        * @return true if rhs and the current functor are equal, false otherwise
        */
        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator==(const ITailOptFunctor<RealType> &rhs) const
        {
            return (this->m_zsize == rhs.m_zsize) &&
                   (this->m_zfirst == rhs.m_zfirst) &&
                   (this->m_zmin == rhs.m_zmin) &&
                   (this->m_zswitch == rhs.m_zswitch) &&
                   (this->m_exp_zswitch == rhs.m_exp_zswitch) &&
                   (this->m_coeff_first == rhs.m_coeff_first) &&
                   (this->m_coeff_second == rhs.m_coeff_second);
        }

        /**
        * Serializes the functor (not suitable for GPUs)
        *
        * @return an std::vector<char> containing the binary representation of the functor
        */
        std::vector<char> serialize() const
        {
            using namespace utils;

            std::vector<char> raw;

            serialization::put_in(m_zsize, raw);
            serialization::put_in(m_zfirst, raw);
            serialization::put_in(m_zmin, raw);
            serialization::put_in(m_zswitch, raw);
            serialization::put_in(m_exp_zswitch, raw);
            serialization::put_in(m_coeff_first, raw);
            serialization::put_in(m_coeff_second, raw);

            return raw;
        }

        /**
        * Deserializes a TailOptFunctor object from binary data (not suitable for GPUs)
        *
        * @tparam CharIter a char interator type
        *
        * @param[in,out] it a char interator
        *
        * @return a LinFunctor<RealType> object
        */
        template <typename CharIter>
        static ITailOptFunctor<RealType> deserialize(CharIter &it)
        {
            using namespace utils;

            auto functor = ITailOptFunctor<RealType>{};

            functor.m_zsize = serialization::get_out<int>(it);
            functor.m_zfirst = serialization::get_out<int>(it);
            functor.m_zmin = serialization::get_out<RealType>(it);
            functor.m_zswitch = serialization::get_out<RealType>(it);
            functor.m_exp_zswitch = serialization::get_out<RealType>(it);
            functor.m_coeff_first = serialization::get_out<RealType>(it);
            functor.m_coeff_second = serialization::get_out<RealType>(it);

            return functor;
        }

    private:
        int m_zsize;
        int m_zfirst;
        RealType m_zmin;
        RealType m_zswitch;
        RealType m_exp_zswitch;
        RealType m_coeff_first;
        RealType m_coeff_second;
    };

}
}
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_DETAIL
