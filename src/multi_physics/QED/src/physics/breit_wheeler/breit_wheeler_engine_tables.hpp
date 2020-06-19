#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES

#include <math.h>
#include <algorithm>
#include <utility>

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Uses picsar tables
#include "../../containers/picsar_tables.hpp"

#include "../../containers/picsar_array.hpp"

#include "../../containers/picsar_span.hpp"

#include "../../math/math_constants.h"

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    enum class dndt_table_type {log};
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
    };

    template<
        typename RealType,
        typename VectorType,
        dndt_table_type TableType = dndt_table_type::log,
        dndt_table_out_policy TableOutPolicy = dndt_table_out_policy::approx
        >
    class dndt_lookup_table{

        public:

            typedef const dndt_lookup_table<
                RealType, containers::picsar_span<RealType>,
                TableType,TableOutPolicy> view_type;

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            dndt_lookup_table(dndt_lookup_table_params<RealType> params);

            std::pair<bool, view_type> get_view();

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType interp(RealType chi_phot) const noexcept;

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            RealType
            interp_flag_out(RealType chi_phot, bool& is_out) const noexcept;

            std::vector<RealType> get_all_coordinates();

            bool set_all_vals(const std::vector<RealType>& vals);

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            bool is_init();
        };

    //Coefficients for che asymptotic behaviour of the dndt table
    //using Erber approximation
    //Erber T. (1966), Reviews of Modern Physics,38, 626
    template <typename T>
    constexpr T erber_dndt_asynt_a = 0.16;
    template <typename T>
    constexpr T erber_dndt_asynt_b = 4./3.;

    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_approx_left(RealType chi_phot)
    {
        using namespace math;
        constexpr auto a = erber_dndt_asynt_a<RealType>;
        constexpr auto b = erber_dndt_asynt_b<RealType>;

        return (pi<RealType>*a/(two<RealType>*b))*chi_phot*chi_phot*
            exp(-two<RealType>*b/chi_phot);
    }

    template <typename RealType>
    PXRMP_INTERNAL_GPU_DECORATOR
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_approx_right(RealType chi_phot)
    {
        using namespace math;
        constexpr auto a = erber_dndt_asynt_a<RealType>;
        constexpr auto b = erber_dndt_asynt_b<RealType>;

        const RealType coeff = tgamma(one<>/three<>)/two<>;
        return a*chi_phot*coeff*
            pow(chi_phot*two<RealType>/b,two_thirds<RealType>);
    }

    template<
        typename RealType,
        typename VectorType,
        dndt_table_out_policy TableOutPolicy
        >
    class dndt_lookup_table<
        RealType, VectorType,
        dndt_table_type::log, TableOutPolicy>{

        public:

            typedef const dndt_lookup_table<
                RealType, containers::picsar_span<RealType>,
                dndt_table_type::log,TableOutPolicy> view_type;

            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            dndt_lookup_table(dndt_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    log(params.chi_phot_min),
                    log(params.chi_phot_max),
                    VectorType(params.chi_phot_how_many)}}
            {};


            PXRMP_INTERNAL_GPU_DECORATOR
            PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
            dndt_lookup_table(dndt_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    log(params.chi_phot_min),
                    log(params.chi_phot_max),
                    vals}}
            {};

            std::pair<bool, view_type> get_view()
            {
                const auto span = containers::picsar_span<RealType>{
                    static_cast<size_t>(m_params.chi_phot_how_many),
                    m_table.m_values.data()
                };
                const view_type view{m_params, span};
                return std::make_pair(m_init_flag, view);
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

        private:
            const dndt_lookup_table_params<RealType> m_params;
            bool m_init_flag = false;
            containers::equispaced_1d_table<RealType, VectorType> m_table;
    };


}
}
}
}

#endif // PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES
