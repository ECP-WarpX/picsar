#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES

//This .hpp file contais the implementation of the lookup tables
// for Breit-Wheeler pair production.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) A.I.Nikishov. & V.I. Ritus Sov. Phys. JETP 19, 2 (1964)
// 2) T.Erber Rev. Mod. Phys. 38, 626 (1966)
// 3) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 4) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar tables
#include "../../containers/picsar_tables.hpp"
//Uses GPU-friendly arrays
#include "../../containers/picsar_array.hpp"
//Uses picsar_span
#include "../../containers/picsar_span.hpp"
//Uses mathematical constants
#include "../../math/math_constants.h"
//Uses serialization
#include "../../utils/serialization.hpp"
//Uses interpolation and upper_bound
#include "../../utils/picsar_algo.hpp"
//Uses cbrt, log and exp
#include "../../math/cmath_overloads.hpp"

#include <algorithm>
#include <vector>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    //Reasonable default values for the dndt_lookup_table_params
    //and the pair_prod_lookup_table_params (see below)
    template <typename T>
    constexpr T default_chi_phot_min = 1.0e-3; /*Default minimum photon chi parameter*/
    template <typename T>
    constexpr T default_chi_phot_max = 1.0e3; /* Default maximum photon chi parameter*/
    const int default_chi_phot_how_many = 256; /* Default number of grid points for photon chi */
    const int default_how_many_frac = 256; /* Default number of grid points for particle chi */

    /**
    * This structure holds the parameters to generate a dN/dt
    * lookup table. The dN/dt lookup table stores the values of the
    * T function (see validation script)
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    struct dndt_lookup_table_params{
        RealType chi_phot_min; /*Minimum photon chi parameter*/
        RealType chi_phot_max;/*Maximum photon chi parameter*/
        int chi_phot_how_many; /* Number of grid points for photon chi */

        /**
        * Operator==
        *
        * @param[in] rhs a structure of the same type
        * @return true if rhs is equal to *this. false otherwise
        */
        bool operator== (const dndt_lookup_table_params<RealType> &rhs) const
        {
            return (chi_phot_min == rhs.chi_phot_min) &&
                (chi_phot_max == rhs.chi_phot_max) &&
                (chi_phot_how_many == rhs.chi_phot_how_many);
        }
    };

    /**
    * The defaultd dndt_lookup_table_params
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    constexpr auto default_dndt_lookup_table_params =
        dndt_lookup_table_params<RealType>{default_chi_phot_min<RealType>,
                                           default_chi_phot_max<RealType>,
                                           default_chi_phot_how_many};
    //__________________________________________________________

    //If a photon has a chi parameter which is out of table,
    //a simple analytical approximation can be used to calculate
    //the T function (see validation script)
    //
    // References:
    // 1) Erber T. (1966), Reviews of Modern Physics,38, 626

    //Coefficients for che asymptotic behaviour of the T function
    template <typename T> //0.16*pi*(3/8)
    constexpr T erber_dndt_asynt_a = 0.1884955592153876; /*Coefficient a of the Erber approximation*/
    template <typename T> //0.16*(gamma(1/3)*(3/2)**(1/3) /2)**2
    constexpr T erber_dndt_asynt_b = 0.37616610710114734; /*Coefficient b of the Erber approximation*/

    /**
    * This function provides an approximation for the T function
    * when chi << 1
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of a photon
    * @return the Erber approximation for T when chi << 1
    */
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_approx_left(RealType chi_phot)
    {
        constexpr RealType coeff = 8./3.;
        return erber_dndt_asynt_a<RealType>*math::m_exp(-coeff/chi_phot);
    }

    /**
    * This function provides an approximation for the T function
    * when chi >> 1
    *
    * @tparam RealType the floating point type to be used
    * @param[in] chi_phot the chi parameter of a photon
    * @return the Erber approximation for T when chi >> 1
    */
    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_approx_right(RealType chi_phot)
    {
        return erber_dndt_asynt_b<RealType>/math::m_cbrt(chi_phot);
    }
    //__________________________________________________________


    /**
    * This class provides the lookup table for dN/dt,
    * storing the values of the T function (see validation script)
    * and providing methods to perform interpolations.
    * It also provides methods for serialization (export to byte array,
    * import from byte array) and to generate "table views" based on
    * non-owning poynters (this is crucial in order to use the table
    * in GPU kernels, as explained below).
    *
    * Internally, this table stores log(T(log(chi))).
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used internally (e.g. std::vector)
    */
    template<
        typename RealType,
        typename VectorType>
    class dndt_lookup_table
    {
        public:

            /**
            * A view_type is essentially a dN/dt lookup table which
            * uses non-owning, constant, pointers to hold the data.
            * The use of views is crucial for GPUs usage. As demonstrated
            * in test_gpu/test_breit_wheeler.cu, it is possible to:
            * - define a thin wrapper around a thrust::device_vector
            * - use this wrapper as a template parameter for a dndt_lookup_table
            * - initialize the lookup table on the CPU
            * - generate a view
            * - pass by copy this view to a GPU kernel (see get_view() method)
            *
            * @tparam RealType the floating point type to be used
            */
            typedef const dndt_lookup_table<
                RealType, containers::picsar_span<const RealType>> view_type;

            /**
            * Constructor (not usable on GPUs).
            * After construction the table is empty. The user has to generate
            * the T function values before being able to use the table.
            *
            * @param params table parameters
            */
            dndt_lookup_table(
                dndt_lookup_table_params<RealType> params = default_dndt_lookup_table_params<RealType>):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    math::m_log(params.chi_phot_min),
                    math::m_log(params.chi_phot_max),
                    VectorType(params.chi_phot_how_many)}}
            {};

            /**
            * Constructor (not usable on GPUs).
            * After construction the table is empty. The user has to generate
            * the T function values before being able to use the table.
            *
            * @param params parameters for table generation
            * @param vals values of the T function
            */
            dndt_lookup_table(dndt_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    math::m_log(params.chi_phot_min),
                    math::m_log(params.chi_phot_max),
                    vals}}
            {
                m_init_flag = true;
            };

            void generate();

            dndt_lookup_table(std::vector<char>& raw_data)
            {
                using namespace utils;

                constexpr size_t min_size =
                    sizeof(char)+//magic_number
                    sizeof(m_params);

                if (raw_data.size() < min_size)
                    throw "Binary data is too small to be a Breit Wheeler \
                     T-function lookup-table.";

                auto it_raw_data = raw_data.begin();

                if (serialization::get_out<char>(it_raw_data) !=
                    static_cast<char>(sizeof(RealType))){
                    throw "Mismatch between RealType used to write and to read \
                        the Breit Wheeler T-function lookup-table";
                }

                m_params = serialization::get_out<
                    dndt_lookup_table_params<RealType>>(it_raw_data);
                m_table = containers::equispaced_1d_table<
                    RealType, VectorType>{std::vector<char>(it_raw_data,
                        raw_data.end())};

                m_init_flag = true;
            };

            PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            bool operator== (
                const dndt_lookup_table<
                    RealType, VectorType> &b) const
            {
                return
                    (m_params == b.m_params) &&
                    (m_init_flag == b.m_init_flag) &&
                    (m_table == b.m_table);
            }

            view_type get_view()
            {
                if(!m_init_flag)
                    throw "Can't generate a view of an uninitialized table";
                const auto span = containers::picsar_span<const RealType>{
                    static_cast<size_t>(m_params.chi_phot_how_many),
                    m_table.get_values_reference().data()
                };
                const view_type view{m_params, span};
                return view;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(RealType chi_phot)  const noexcept
            {
                if (chi_phot < m_params.chi_phot_min){
                    return dndt_approx_left<RealType>(chi_phot);
                }

                if(chi_phot > m_params.chi_phot_max){
                    return dndt_approx_right<RealType>(chi_phot);
                }
                return math::m_exp(m_table.interp(math::m_log(chi_phot)));
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp_flag_out(
                const RealType chi_phot, bool& is_out) const noexcept
            {
                is_out = false;
                if(chi_phot < m_params.chi_phot_min || chi_phot >  m_params.chi_phot_max)
                    is_out = true;
                return interp(chi_phot);
            }

            std::vector<RealType> get_all_coordinates() const noexcept
            {
                auto all_coords = m_table.get_all_coordinates();
                std::transform(all_coords.begin(),all_coords.end(),all_coords.begin(),
                    [](RealType a){return math::m_exp(a);});
                return all_coords;
            }

            bool set_all_vals(const std::vector<RealType>& vals)
            {
                if(vals.size() == m_table.get_how_many_x()){
                    for(int i = 0; i < vals.size(); ++i){
                        m_table.set_val(i, math::m_log(vals[i]));
                    }
                    m_init_flag = true;
                    return true;
                }
                return false;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            bool is_init()
            {
                return m_init_flag;
            }

            std::vector<char> serialize()
            {
                using namespace utils;

                if(!m_init_flag)
                    throw "Cannot serialize an unitialized table";

                std::vector<char> res;

                serialization::put_in(static_cast<char>(sizeof(RealType)), res);
                serialization::put_in(m_params, res);

                auto tdata = m_table.serialize();
                res.insert(res.end(), tdata.begin(), tdata.end());

                return res;
            }

        private:
            dndt_lookup_table_params<RealType> m_params;
            bool m_init_flag = false;
            containers::equispaced_1d_table<RealType, VectorType> m_table;
    };

    //________________________________________________________________________________

    template<typename RealType>
    struct pair_prod_lookup_table_params{
        RealType chi_phot_min; //Min chi_phot
        RealType chi_phot_max; //Max chi_phot
        int chi_phot_how_many; //How many points
        int how_many_frac; //How many points

        bool operator== (
            const pair_prod_lookup_table_params<RealType> &b) const
        {
            return (chi_phot_min == b.chi_phot_min) &&
                (chi_phot_max == b.chi_phot_max) &&
                (chi_phot_how_many == b.chi_phot_how_many) &&
                (how_many_frac == b.how_many_frac);
        }
    };

        template<
        typename RealType,
        typename VectorType>
    class pair_prod_lookup_table_logchi_linfrac{
        public:
            typedef const pair_prod_lookup_table_logchi_linfrac<
                RealType, containers::picsar_span<const RealType>> view_type;

            pair_prod_lookup_table_logchi_linfrac(
                pair_prod_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_2d_table<RealType, VectorType>{
                    math::m_log(params.chi_phot_min),
                    math::m_log(params.chi_phot_max),
                    math::zero<RealType>,
                    math::half<RealType>,
                    params.chi_phot_how_many, params.how_many_frac,
                    VectorType(params.chi_phot_how_many * params.how_many_frac)}}
            {};

            pair_prod_lookup_table_logchi_linfrac(
                pair_prod_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_2d_table<RealType, VectorType>{
                    math::m_log(params.chi_phot_min),
                    math::m_log(params.chi_phot_max),
                    math::zero<RealType>,
                    math::half<RealType>,
                    params.chi_phot_how_many, params.how_many_frac,
                    vals}}
            {
                m_init_flag = true;
            };

            pair_prod_lookup_table_logchi_linfrac(std::vector<char>& raw_data)
            {
                using namespace utils;

                constexpr size_t min_size =
                    sizeof(char)+//magic_number
                    sizeof(m_params);

                if (raw_data.size() < min_size)
                    throw "Binary data is too small to be a Breit Wheeler \
                     pair production lookup-table.";

                auto it_raw_data = raw_data.begin();

                if (serialization::get_out<char>(it_raw_data) !=
                    static_cast<char>(sizeof(RealType))){
                    throw "Mismatch between RealType used to write and to read \
                        the Breit pair production lookup-table";
                }

                m_params = serialization::get_out<
                    pair_prod_lookup_table_params<RealType>>(it_raw_data);
                m_table = containers::equispaced_2d_table<
                    RealType, VectorType>{std::vector<char>(it_raw_data,
                        raw_data.end())};

                m_init_flag = true;
            };

            PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            bool operator== (
                const pair_prod_lookup_table_logchi_linfrac<
                    RealType, VectorType> &b) const
            {
                return
                    (m_params == b.m_params) &&
                    (m_init_flag == b.m_init_flag) &&
                    (m_table == b.m_table);
            }

            view_type get_view()
            {
                if(!m_init_flag)
                    throw "Can't generate a view of an uninitialized table";
                const auto span = containers::picsar_span<const RealType>{
                    static_cast<size_t>(m_params.chi_phot_how_many *
                        m_params.how_many_frac),
                        m_table.get_values_reference().data()
                };
                const view_type view{m_params, span};
                return view;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(
                const RealType chi_phot,
                const RealType unf_zero_one_minus_epsi) const noexcept
            {
                using namespace math;

                auto e_chi_phot = chi_phot;
                if(chi_phot<m_params.chi_phot_min)
                    e_chi_phot = m_params.chi_phot_min;
                else if (chi_phot > m_params.chi_phot_max)
                    e_chi_phot = m_params.chi_phot_max;
                const auto log_e_chi_phot = m_log(e_chi_phot);

                const auto prob = unf_zero_one_minus_epsi*half<RealType>;

                const auto upper_frac_index = utils::picsar_upper_bound_functor(
                    0, m_params.how_many_frac,prob,[&](int i){
                        return (m_table.interp_first_coord(
                            log_e_chi_phot, i));
                    });
                const auto lower_frac_index = upper_frac_index-1;

                const auto upper_frac = m_table.get_y_coord(upper_frac_index);
                const auto lower_frac = m_table.get_y_coord(lower_frac_index);

                const auto lower_prob= m_table.interp_first_coord
                    (log_e_chi_phot, lower_frac_index);
                const auto upper_prob = m_table.interp_first_coord
                    (log_e_chi_phot, upper_frac_index);

                const auto frac = utils::linear_interp(
                    lower_prob, upper_prob, lower_frac, upper_frac,
                    prob);

                const auto chi = frac*chi_phot;

                auto res = (unf_zero_one_minus_epsi < half<RealType>)?chi:(chi_phot - chi);

                return (unf_zero_one_minus_epsi < half<RealType>)?
                    chi:(chi_phot - chi);
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp_flag_out(
                const RealType chi_phot,
                const RealType unf_zero_one_minus_epsi,
                bool& is_out) const noexcept
            {
                is_out = false;
                if(chi_phot < m_params.chi_phot_min || chi_phot >  m_params.chi_phot_max)
                    is_out = true;
                return interp(chi_phot, unf_zero_one_minus_epsi);
            }

            std::vector<std::array<RealType,2>> get_all_coordinates() const noexcept
            {
                auto all_coords = m_table.get_all_coordinates();
                std::transform(all_coords.begin(),all_coords.end(),all_coords.begin(),
                    [](std::array<RealType,2> a){return
                        std::array<RealType,2>{math::m_exp(a[0]), a[1]};});
                return all_coords;
            }

            bool set_all_vals(const std::vector<RealType>& vals)
            {
                if(vals.size() == m_table.get_how_many_x()*
                    m_table.get_how_many_y()){
                    for(int i = 0; i < vals.size(); ++i){
                        m_table.set_val(i,vals[i]);
                    }
                    m_init_flag = true;
                    return true;
                }
                return false;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            bool is_init()
            {
                return m_init_flag;
            }

            std::vector<char> serialize()
            {
                using namespace utils;

                if(!m_init_flag)
                    throw "Cannot serialize an unitialized table";

                std::vector<char> res;

                serialization::put_in(static_cast<char>(sizeof(RealType)), res);
                serialization::put_in(m_params, res);

                auto tdata = m_table.serialize();
                res.insert(res.end(), tdata.begin(), tdata.end());

                return res;
            }

        private:
            pair_prod_lookup_table_params<RealType> m_params;
            bool m_init_flag = false;
            containers::equispaced_2d_table<RealType, VectorType> m_table;
    };

}
}
}
}

#endif // PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
