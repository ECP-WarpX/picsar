#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES

#include <math.h>
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

#include "../../utils/picsar_algo.hpp"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    enum class dndt_table_out_policy {approx, extrema};

    template<typename RealType>
    struct dndt_lookup_table_params{
        // dN/dt table:
        //__breit_wheeler_min_tdndt_chi_phot  is the inferior
        //limit of the total pair production rate lookup table.
        //If  __breit_wheeler_min_chi_phot < chi <__breit_wheeler_min_tdndt_chi_phot
        //BW process is taken into account, but the Erber approximation
        //rather than the lookup table is used.
        RealType chi_phot_min; //Min chi_phot
        RealType chi_phot_max; //Max chi_phot
        int chi_phot_how_many; //How many points

        bool operator== (const dndt_lookup_table_params<RealType> &b) const
        {
            return (chi_phot_min == b.chi_phot_min) &&
                (chi_phot_max == b.chi_phot_max) &&
                (chi_phot_how_many == b.chi_phot_how_many);
        }
    };


    //Coefficients for che asymptotic behaviour of the dndt table
    //using Erber approximation
    //Erber T. (1966), Reviews of Modern Physics,38, 626
    template <typename T> //0.16*pi*(3/8)
    constexpr T erber_dndt_asynt_a = 0.1884955592153876;
    template <typename T> //0.16*(gamma(1/3)*(3/2)**(1/3) /2)**2
    constexpr T erber_dndt_asynt_b = 0.37616610710114734;


    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_approx_left(RealType chi_phot)
    {
        constexpr RealType coeff = 8./3.;
        return erber_dndt_asynt_a<RealType>*exp(-coeff/chi_phot);
    }

    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_approx_right(RealType chi_phot)
    {
        return erber_dndt_asynt_b<RealType>/cbrt(chi_phot);
    }

    template<
        typename RealType,
        typename VectorType,
        dndt_table_out_policy TableOutPolicy = dndt_table_out_policy::approx
        >
    class dndt_lookup_table{

        public:

            typedef const dndt_lookup_table<
                RealType, containers::picsar_span<RealType>,
                TableOutPolicy> view_type;

            dndt_lookup_table(dndt_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    log(params.chi_phot_min),
                    log(params.chi_phot_max),
                    VectorType(params.chi_phot_how_many)}}
            {};

            dndt_lookup_table(dndt_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    log(params.chi_phot_min),
                    log(params.chi_phot_max),
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
                const auto span = containers::picsar_span<RealType>{
                    static_cast<size_t>(m_params.chi_phot_how_many),
                    m_table.m_values.data()
                };
                const view_type view{m_params, span};
                return view;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(RealType chi_phot)  const noexcept
            {
                PXRMP_INTERNAL_CONSTEXPR_IF
                    (TableOutPolicy == dndt_table_out_policy::extrema){
                        if(chi_phot<m_params.chi_phot_min)
                            chi_phot = m_params.chi_phot_min;
                        else if (chi_phot > m_params.chi_phot_max)
                            chi_phot = m_params.chi_phot_max;
                }
                else
                {
                    if (chi_phot < m_params.chi_phot_min){
                        return dndt_approx_left<RealType>(chi_phot);
                    }

                    if(chi_phot > m_params.chi_phot_max){
                        return dndt_approx_right<RealType>(chi_phot);
                    }
                }

                return exp(m_table.interp(log(chi_phot)));
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
                    [](RealType a){return exp(a);});
                return all_coords;
            }

            bool set_all_vals(const std::vector<RealType>& vals)
            {
                if(vals.size() == m_table.get_how_many_x()){
                    for(int i = 0; i < vals.size(); ++i){
                        m_table.set_val(i, log(vals[i]));
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
                RealType, containers::picsar_span<RealType>> view_type;

            pair_prod_lookup_table_logchi_linfrac(
                pair_prod_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_2d_table<RealType, VectorType>{
                    log(params.chi_phot_min),
                    log(params.chi_phot_max),
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
                    log(params.chi_phot_min),
                    log(params.chi_phot_max),
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
                const auto span = containers::picsar_span<RealType>{
                    static_cast<size_t>(m_params.chi_phot_how_many),
                    m_table.m_values.data()
                };
                const view_type view{m_params, span};
                return view;
            }

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(
                RealType chi_phot,
                const RealType unf_zero_one_minus_epsi) const noexcept
            {
                using namespace math;

                if(chi_phot<m_params.chi_phot_min)
                    chi_phot = m_params.chi_phot_min;
                else if (chi_phot > m_params.chi_phot_max)
                    chi_phot = m_params.chi_phot_max;

                const auto prob = unf_zero_one_minus_epsi*half<RealType>;

                const auto upper_frac_index = utils::picsar_upper_bound_functor(
                    0, m_params.how_many_frac,prob,[&](int i){
                        return (m_table.interp_first_coord(
                            log(chi_phot), i));
                    });
                const auto lower_frac_index = upper_frac_index-1;

                const auto upper_frac = m_table.get_y_coord(upper_frac_index);
                const auto lower_frac = m_table.get_y_coord(lower_frac_index);

                const auto lower_prob= m_table.interp_first_coord
                    (log(chi_phot), lower_frac_index);
                const auto upper_prob = m_table.interp_first_coord
                    (log(chi_phot), upper_frac_index);

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
                        std::array<RealType,2>{exp(a[0]), a[1]};});
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
