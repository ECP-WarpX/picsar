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
        LinFunctor(){};

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
        * @tparam Iter a char interator
        * @return a LinFunctor<RealType> object
        */
        template <typename Iter>
        static LinFunctor<RealType> deserialize(Iter &it)
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
    * The functor is essentially the inverse of LinFunctor.
    * It also contains methods to serialize and deserialize itself.
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
        ILinFunctor(){};

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
        * @param[in] rhs a const reference to a ILinFunctor<RealType>
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
        * @tparam Iter a char interator
        * @return a ILinFunctor<RealType> object
        */
        template <typename Iter>
        static ILinFunctor<RealType> deserialize(Iter &it)
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

    template <typename RealType>
    class TailOptFunctor
    {
    public:

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        TailOptFunctor(){};

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        TailOptFunctor(const int zsize, const int zfirst,
                       const RealType zmin, const RealType zmax, const RealType zswitch) : m_zsize{zsize}, m_zfirst{zfirst},
                                                                                           m_zmin{zmin}, m_zmax{zmax}, m_zswitch{zswitch}
        {
            m_logzmin = math::m_log(m_zmin);
            m_logzswitch = math::m_log(m_zswitch);
        }

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        RealType operator()(const int i) const
        {
            using namespace picsar::multi_physics::math;

            if (i < m_zfirst)
            {
                return m_logzmin + i * (m_logzswitch - m_logzmin) / (m_zfirst - 1);
            }
            else
            {
                return math::m_log(
                    m_zswitch + ((i + 1) - m_zfirst) * (m_zmax - m_zswitch) / (m_zsize - m_zfirst));
            }
        }

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator==(const TailOptFunctor<RealType> &rhs) const
        {
            return (this->m_zsize == rhs.m_zsize) &&
                   (this->m_zfirst == rhs.m_zfirst) &&
                   (this->m_zmin == rhs.m_zmin) &&
                   (this->m_logzmin == rhs.m_logzmin) &&
                   (this->m_zmax == rhs.m_zmax) &&
                   (this->m_zswitch == rhs.m_zswitch) &&
                   (this->m_logzswitch == rhs.m_logzswitch);
        }

        std::vector<char> serialize() const
        {
            using namespace utils;

            std::vector<char> raw;

            serialization::put_in(m_zsize, raw);
            serialization::put_in(m_zfirst, raw);
            serialization::put_in(m_zmin, raw);
            serialization::put_in(m_logzmin, raw);
            serialization::put_in(m_zmax, raw);
            serialization::put_in(m_zswitch, raw);
            serialization::put_in(m_logzswitch, raw);

            return raw;
        }

        template <typename Iter>
        static TailOptFunctor<RealType> deserialize(Iter &it)
        {
            using namespace utils;

            auto functor = TailOptFunctor<RealType>{};

            functor.m_zsize = serialization::get_out<int>(it);
            functor.m_zfirst = serialization::get_out<int>(it);
            functor.m_zmin = serialization::get_out<RealType>(it);
            functor.m_logzmin = serialization::get_out<RealType>(it);
            functor.m_zmax = serialization::get_out<RealType>(it);
            functor.m_zswitch = serialization::get_out<RealType>(it);
            functor.m_logzswitch = serialization::get_out<RealType>(it);

            return functor;
        }

    private:
        int m_zsize;
        int m_zfirst;
        RealType m_zmin;
        RealType m_logzmin;
        RealType m_zmax;
        RealType m_zswitch;
        RealType m_logzswitch;
    };

    template <typename RealType>
    class ITailOptFunctor
    {
    public:

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        ITailOptFunctor(){};

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        ITailOptFunctor(const int zsize, const int zfirst,
                        const RealType zmin, const RealType zmax, const RealType zswitch) : m_zsize{zsize}, m_zfirst{zfirst},
                                                                                            m_zmin{zmin}, m_zmax{zmax}, m_zswitch{zswitch}
        {
            m_logzmin = math::m_log(m_zmin);
            m_logzswitch = math::m_log(m_zswitch);
        }

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        int operator()(const RealType logz) const
        {
            using namespace picsar::multi_physics::math;

            if (logz < m_logzswitch)
            {
                return static_cast<int>(
                    math::m_round((m_zfirst - 1) * (logz - m_logzmin) / (m_logzswitch - m_logzmin)));
            }
            else
            {
                const auto z = m_exp(logz);
                return static_cast<int>(
                    m_floor(m_zfirst + (m_zsize - m_zfirst - 1) * (z - m_zswitch) / (m_zmax - m_zswitch)));
            }
        }

        PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
        bool operator==(const ITailOptFunctor<RealType> &rhs) const
        {
            return (this->m_zsize == rhs.m_zsize) &&
                   (this->m_zfirst == rhs.m_zfirst) &&
                   (this->m_zmin == rhs.m_zmin) &&
                   (this->m_logzmin == rhs.m_logzmin) &&
                   (this->m_zmax == rhs.m_zmax) &&
                   (this->m_zswitch == rhs.m_zswitch) &&
                   (this->m_logzswitch == rhs.m_logzswitch);
        }

        std::vector<char> serialize() const
        {
            using namespace utils;

            std::vector<char> raw;

            serialization::put_in(m_zsize, raw);
            serialization::put_in(m_zfirst, raw);
            serialization::put_in(m_zmin, raw);
            serialization::put_in(m_logzmin, raw);
            serialization::put_in(m_zmax, raw);
            serialization::put_in(m_zswitch, raw);
            serialization::put_in(m_logzswitch, raw);

            return raw;
        }

        template <typename Iter>
        static ITailOptFunctor<RealType> deserialize(Iter &it)
        {
            using namespace utils;

            auto functor = ITailOptFunctor<RealType>{};

            functor.m_zsize = serialization::get_out<int>(it);
            functor.m_zfirst = serialization::get_out<int>(it);
            functor.m_zmin = serialization::get_out<RealType>(it);
            functor.m_logzmin = serialization::get_out<RealType>(it);
            functor.m_zmax = serialization::get_out<RealType>(it);
            functor.m_zswitch = serialization::get_out<RealType>(it);
            functor.m_logzswitch = serialization::get_out<RealType>(it);

            return functor;
        }

    private:
        int m_zsize;
        int m_zfirst;
        RealType m_zmin;
        RealType m_logzmin;
        RealType m_zmax;
        RealType m_zswitch;
        RealType m_logzswitch;
    };
}
}
}
}
}

#endif //PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_DETAIL
