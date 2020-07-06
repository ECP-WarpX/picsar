#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES

#include <algorithm>
#include <utility>
#include <vector>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar tables
#include "../../containers/picsar_tables.hpp"
#include "../../containers/picsar_array.hpp"
#include "../../containers/picsar_span.hpp"

#include "../../math/math_constants.h"

#include "../../utils/serialization.hpp"
//Uses exp and log
#include "../../math/cmath_overloads.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    enum class dndt_table_out_policy {approx, extrema};

    template<typename RealType>
    struct dndt_lookup_table_params{
        RealType chi_part_min; //Min chi_phot
        RealType chi_part_max; //Max chi_phot
        int chi_part_how_many; //How many points

        bool operator== (const dndt_lookup_table_params<RealType> &b) const
        {
            return (chi_part_min == b.chi_part_min) &&
                (chi_part_max == b.chi_part_max) &&
                (chi_part_how_many == b.chi_part_how_many);
        }
    };

    template<
        typename RealType,
        typename VectorType,
        dndt_table_out_policy TableOutPolicy
        >
    class dndt_lookup_table{

        public:

            typedef const dndt_lookup_table<
                RealType, containers::picsar_span<const RealType>,
                TableOutPolicy> view_type;

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
                const dndt_lookup_table<
                    RealType, VectorType, TableOutPolicy> &b) const
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
                PXRMP_INTERNAL_CONSTEXPR_IF
                    (TableOutPolicy == dndt_table_out_policy::extrema){
                        if(chi_part<m_params.chi_part_min)
                            chi_part = m_params.chi_part_min;
                        else if (chi_part > m_params.chi_part_max)
                            chi_part = m_params.chi_part_max;
                }
                else
                {
                    if (chi_part < m_params.chi_part_min){
                        /*TODO*/
                        return math::zero<RealType>;
                    }

                    if(chi_part > m_params.chi_part_max){
                        /*TODO*/
                        return  math::zero<RealType>;
                    }
                }

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

        private:
            dndt_lookup_table_params<RealType> m_params;
            bool m_init_flag = false;
            containers::equispaced_1d_table<RealType, VectorType> m_table;
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
class photon_emission_lookup_table_logchi_logfrac
{
public:

    typedef const photon_emission_lookup_table_logchi_logfrac<
        RealType, containers::picsar_span<const RealType>> view_type;

    photon_emission_lookup_table_logchi_logfrac(
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


    photon_emission_lookup_table_logchi_logfrac(
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

    photon_emission_lookup_table_logchi_logfrac(std::vector<char>& raw_data)
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
        const photon_emission_lookup_table_logchi_logfrac<
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

private:
    photon_emission_lookup_table_params<RealType> m_params;
    bool m_init_flag = false;
    containers::equispaced_2d_table<RealType, VectorType> m_table;
};


}
}
}
}

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
