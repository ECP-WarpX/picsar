#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES

//This .hpp file contais the implementation of the lookup tables
//for Quantum Synchrotron photon production.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

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
//Uses log and exp
#include "../../math/cmath_overloads.hpp"

#include <algorithm>
#include <vector>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    //Reasonable default values for the dndt_lookup_table_params
    //and the pair_prod_lookup_table_params (see below)
    template <typename T>
    constexpr T default_chi_part_min = 1.0e-3; /*Default minimum particle chi parameter*/
    template <typename T>
    constexpr T default_chi_part_max = 1.0e3; /* Default maximum particle chi parameter*/
    const int default_chi_part_how_many = 256; /* Default number of grid points for particle chi */
    const int default_how_many_frac = 256; /* Default number of grid points for photon chi */

    /**
    * This structure holds the parameters to generate a dN/dt
    * lookup table. The dN/dt lookup table stores the values of the
    * T function (see validation script)
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    struct dndt_lookup_table_params{
        RealType chi_part_min; /*Minimum particle chi parameter*/
        RealType chi_part_max; /*Maximum particle chi parameter*/
        int chi_part_how_many; /* Number of grid points for particle chi */

        /**
        * Operator==
        *
        * @param[in] rhs a structure of the same type
        * @return true if rhs is equal to *this. false otherwise
        */
        bool operator== (const dndt_lookup_table_params<RealType> &rhs) const
        {
            return (chi_part_min == rhs.chi_part_min) &&
                (chi_part_max == rhs.chi_part_max) &&
                (chi_part_how_many == rhs.chi_part_how_many);
        }
    };

    /**
    * The default dndt_lookup_table_params
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    constexpr auto default_dndt_lookup_table_params =
        dndt_lookup_table_params<RealType>{default_chi_part_min<RealType>,
                                           default_chi_part_max<RealType>,
                                           default_chi_part_how_many};
    //__________________________________________________________

    /**
    * generation_policy::force_internal_double can be used to force the
    * calulcations of a lookup tables using double precision, even if
    * the final result is stored in a single precision variable.
    */
    enum class generation_policy{
        regular,
        force_internal_double
    };


    template<
        typename RealType,
        typename VectorType>
    class dndt_lookup_table{

        public:

            typedef const dndt_lookup_table<
                RealType, containers::picsar_span<const RealType>> view_type;

            dndt_lookup_table(dndt_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    math::m_log(params.chi_part_min),
                    math::m_log(params.chi_part_max),
                    VectorType(params.chi_part_how_many)}}
            {};

            dndt_lookup_table(dndt_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    math::m_log(params.chi_part_min),
                    math::m_log(params.chi_part_max),
                    vals}}
            {
                m_init_flag = true;
            };

            template <generation_policy Policy = generation_policy::regular>
            void generate(const bool show_progress  = true);

            dndt_lookup_table(std::vector<char>& raw_data)
            {
                using namespace utils;

                constexpr size_t min_size =
                    sizeof(char)+//magic_number
                    sizeof(m_params);

                if (raw_data.size() < min_size)
                    throw "Binary data is too small to be a Quantum Synchrotron \
                     G-function lookup-table.";

                auto it_raw_data = raw_data.begin();

                if (serialization::get_out<char>(it_raw_data) !=
                    static_cast<char>(sizeof(RealType))){
                    throw "Mismatch between RealType used to write and to read \
                        the Quantum Synchrotron G-function lookup-table";
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
                const dndt_lookup_table<RealType, VectorType> &b) const
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
                    static_cast<size_t>(m_params.chi_part_how_many),
                    m_table.get_values_reference().data()
                };
                const view_type view{m_params, span};
                return view;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(RealType chi_part)  const noexcept
            {
                if(chi_part<m_params.chi_part_min)
                    chi_part = m_params.chi_part_min;
                else if (chi_part > m_params.chi_part_max)
                    chi_part = m_params.chi_part_max;

                return math::m_exp(m_table.interp(math::m_log(chi_part)));
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp_flag_out(
                const RealType chi_part, bool& is_out) const noexcept
            {
                is_out = false;
                if(chi_part < m_params.chi_part_min || chi_part >  m_params.chi_part_max)
                    is_out = true;
                return interp(chi_part);
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

        protected:
            dndt_lookup_table_params<RealType> m_params;
            bool m_init_flag = false;
            containers::equispaced_1d_table<RealType, VectorType> m_table;

        private:
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            static RealType aux_generate_double(RealType x);

    };

    //________________________________________________________________________________

    template<typename RealType>
    struct photon_emission_lookup_table_params{
        RealType chi_part_min; //Min chi_part
        RealType chi_part_max; //Max chi_part
        RealType frac_min; //Max chi_part
        int chi_part_how_many; //How many points
        int how_many_frac; //How many points

        bool operator== (
            const photon_emission_lookup_table_params<RealType> &b) const
        {
            return (chi_part_min == b.chi_part_min) &&
                (chi_part_max == b.chi_part_max) &&
                (frac_min == b.frac_min) &&
                (chi_part_how_many == b.chi_part_how_many) &&
                (how_many_frac == b.how_many_frac);
        }
    };

template<typename RealType, typename VectorType>
class photon_emission_lookup_table
{
public:

    typedef const photon_emission_lookup_table<
        RealType, containers::picsar_span<const RealType>> view_type;

    photon_emission_lookup_table(
            photon_emission_lookup_table_params<RealType> params):
        m_params{params},
        m_table{containers::equispaced_2d_table<RealType, VectorType>{
                math::m_log(params.chi_part_min),
                math::m_log(params.chi_part_max),
                math::m_log(params.frac_min),
                math::m_log(math::one<RealType>),
                params.chi_part_how_many, params.how_many_frac,
                VectorType(params.chi_part_how_many * params.how_many_frac)}}
        {};


    photon_emission_lookup_table(
        photon_emission_lookup_table_params<RealType> params,
        VectorType vals):
    m_params{params},
    m_table{containers::equispaced_2d_table<RealType, VectorType>{
            math::m_log(params.chi_part_min),
            math::m_log(params.chi_part_max),
            math::m_log(params.frac_min),
            math::m_log(math::one<RealType>),
            params.chi_part_how_many, params.how_many_frac,
            vals}}
    {
        m_init_flag = true;
    };

    template <generation_policy Policy = generation_policy::regular>
    void generate(bool show_progress = true);

    photon_emission_lookup_table(std::vector<char>& raw_data)
    {
        using namespace utils;

        constexpr size_t min_size =
            sizeof(char)+//magic_number
            sizeof(m_params);

        if (raw_data.size() < min_size)
            throw "Binary data is too small to be a Quantum Synchrotron \
            emisson lookup-table.";

        auto it_raw_data = raw_data.begin();

        if (serialization::get_out<char>(it_raw_data) !=
            static_cast<char>(sizeof(RealType))){
            throw "Mismatch between RealType used to write and to read \
                the Quantum Synchrotron lookup-table";
        }

        m_params = serialization::get_out<
            photon_emission_lookup_table_params<RealType>>(it_raw_data);
        m_table = containers::equispaced_2d_table<
            RealType, VectorType>{std::vector<char>(it_raw_data,
                raw_data.end())};

        m_init_flag = true;
    };

    PXRMP_INTERNAL_GPU_DECORATOR PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    bool operator== (
        const photon_emission_lookup_table<
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
            static_cast<size_t>(m_params.chi_part_how_many *
                m_params.how_many_frac), m_table.get_values_reference().data()
        };

        const view_type view{m_params, span};
        return view;
    }

    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType interp(
        const RealType chi_part,
        const RealType unf_zero_one_minus_epsi) const noexcept
    {
        using namespace math;

        auto e_chi_part = chi_part;
        if(chi_part<m_params.chi_part_min)
            e_chi_part = m_params.chi_part_min;
        else if (chi_part > m_params.chi_part_max)
            e_chi_part = m_params.chi_part_max;
        const auto log_e_chi_part = m_log(e_chi_part);

        const auto log_prob = m_log(one<RealType>-unf_zero_one_minus_epsi);

        const auto upper_frac_index = utils::picsar_upper_bound_functor(
            0, m_params.how_many_frac,log_prob,[&](int i){
                return (m_table.interp_first_coord(
                    log_e_chi_part, i));
            });

        if(upper_frac_index == 0)
            return m_params.frac_min*chi_part;

        if(upper_frac_index ==  m_params.how_many_frac)
            return chi_part;

        const auto lower_frac_index = upper_frac_index-1;

        const auto upper_log_frac = m_table.get_y_coord(upper_frac_index);
        const auto lower_log_frac = m_table.get_y_coord(lower_frac_index);

        const auto lower_log_prob= m_table.interp_first_coord
            (log_e_chi_part, lower_frac_index);
        const auto upper_log_prob = m_table.interp_first_coord
            (log_e_chi_part, upper_frac_index);

        const auto log_frac = utils::linear_interp(
            lower_log_prob, upper_log_prob, lower_log_frac, upper_log_frac,
            log_prob);
        return  m_exp(log_frac)*chi_part;
    }

    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType interp_flag_out(
        const RealType chi_part,
        const RealType unf_zero_one_minus_epsi,
        bool& is_out) const noexcept
    {
        is_out = false;
        if(chi_part < m_params.chi_part_min || chi_part >  m_params.chi_part_max)
            is_out = true;
        return interp(chi_part, unf_zero_one_minus_epsi);
    }

    std::vector<std::array<RealType,2>> get_all_coordinates() const noexcept
    {
        auto all_coords = m_table.get_all_coordinates();
        std::transform(all_coords.begin(),all_coords.end(),all_coords.begin(),
            [](std::array<RealType,2> a){return
                std::array<RealType,2>{math::m_exp(a[0]), math::m_exp(a[1])};});
        return all_coords;
    }

    bool set_all_vals(const std::vector<RealType>& vals)
    {
        if(vals.size() == m_table.get_how_many_x()*
            m_table.get_how_many_y()){
            for(int i = 0; i < vals.size(); ++i){
                m_table.set_val(i,math::m_log(vals[i]));
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

protected:
    photon_emission_lookup_table_params<RealType> m_params;
    bool m_init_flag = false;
    containers::equispaced_2d_table<RealType, VectorType> m_table;

private:
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    static std::vector<RealType>
    aux_generate_double(RealType x,
        const std::vector<RealType>& y);
};


}
}
}
}

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
