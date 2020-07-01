#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES

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
                RealType, containers::picsar_span<RealType>,
                TableOutPolicy> view_type;

            dndt_lookup_table(dndt_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    log(params.chi_part_min),
                    log(params.chi_part_max),
                    VectorType(params.chi_part_how_many)}}
            {};

            dndt_lookup_table(dndt_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    log(params.chi_part_min),
                    log(params.chi_part_max),
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
                const auto span = containers::picsar_span<RealType>{
                    static_cast<size_t>(m_params.chi_part_how_many),
                    m_table.m_values.data()
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

                return exp(m_table.interp(log(chi_part)));
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


}
}
}
}

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
