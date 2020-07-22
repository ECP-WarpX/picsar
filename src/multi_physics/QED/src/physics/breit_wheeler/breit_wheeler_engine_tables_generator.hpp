#ifndef PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES_GENERATOR
#define PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES_GENERATOR

//This .hpp file extends the implementation of the lookup tables
//for Breit-Wheeler pair production with methods to generate
//the lookup tables.
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

//Implements methods of BW lookup tables
#include "breit_wheeler_engine_tables.hpp"
//Uses BW tabulated functions
#include "breit_wheeler_engine_tabulated_functions.hpp"
//Uses cmath overloads
#include "../../math/cmath_overloads.hpp"
//Uses progress bar
#include "../../utils/progress_bar.hpp"

#include <omp.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <type_traits>
#include <stdexcept>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace breit_wheeler{

    //________________ dN/dt table _____________________________________________

    /**
    * Auxiliary function used to compute the T function in double precision
    * with a single precision argument and to cast back the result to single precision
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used (relevant for the class of which is a method is member)
    *
    * @param[in] x the value to be passed to compute_T_function
    *
    * @return the result of compute_T_function
    */
    template<typename RealType, typename VectorType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_lookup_table<RealType, VectorType>::
    aux_generate_double(RealType x)
    {
        return static_cast<RealType>(
            compute_T_function<double>(x));
    }

    /**
    * Generates the lookup table (not usable on GPUs).
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used (relevant for the class of which is a method is member)
    * @tparam Policy if set to generation_policy::force_internal_double it forces internal calculations in double precision
    *
    * @param[in] show_progress if true it shows a nice progress bar
    */
    template<typename RealType, typename VectorType>
    template<generation_policy Policy>
    void dndt_lookup_table<RealType, VectorType>::generate(
        const bool show_progress)
    {
        constexpr bool use_internal_double =
            (Policy == generation_policy::force_internal_double) &&
            !std::is_same<RealType,double>();

        auto t_start =  std::chrono::system_clock::now();

        const auto all_coords = get_all_coordinates();
        auto all_vals = std::vector<RealType>(all_coords.size());

        int count = 0;
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(all_coords.size()); ++i){
            PXRMP_INTERNAL_CONSTEXPR_IF (use_internal_double){
                all_vals[i] = aux_generate_double(all_coords[i]);
            }
            else {
                all_vals[i] = compute_T_function(all_coords[i]);
            }

            #pragma omp critical
            {
                count++;
                utils::draw_progress(count,
                    all_vals.size(), "Breit-Wheeler dN/dt", 1);
            }
        }

        for (auto& val : all_vals){
            if(std::isnan(val))
                throw std::runtime_error("Error: nan detected in generated table!");
        }

        set_all_vals(all_vals);

        auto t_end =  std::chrono::system_clock::now();
        utils::draw_progress(
            count, all_vals.size(), "Breit-Wheeler dN/dt", 1, true);
        std::cout << " Done in " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                t_end - t_start).count()/1000.0 << " seconds. \n" << std::endl;

        m_init_flag = true;
    }

    //__________________________________________________________________________

    //________________ Pair production table ___________________________________

    /**
    * Auxiliary function used to compute cumulative probability distribution in double precision
    * with a single precision argument and to cast back the result to single precision
    *
    * @tparam RealType the floating point type to be used
    * @tparam RealType the vector type to be used (relevant for the class of which is a method is member)
    *
    * @param[in] x the value to be passed to compute_cumulative_prob
    *
    * @return the result of compute_cumulative_prob
    */
    template<typename RealType, typename VectorType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    std::vector<RealType>
    pair_prod_lookup_table<RealType, VectorType>::
    aux_generate_double(const RealType x, const std::vector<RealType>& y)
    {
        auto dtemp = std::vector<double>(y.size());
        std::transform(y.begin(), y.end(), dtemp.begin(),
            [](RealType yy){ return static_cast<double>(yy); });
        const auto dres = compute_cumulative_prob_opt<
            double, std::vector<double>>(x,dtemp);
        auto res = std::vector<RealType>(y.size());
        std::transform(dres.begin(), dres.end(), res.begin(),
            [](double rr){ return static_cast<RealType>(rr); });
        return res;
    }

    /**
    * Generates the content of the lookup table (not usable on GPUs).
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used (relevant for the class of which is a method is member)
    * @tparam Policy if set to generation_policy::force_internal_double it forces internal calculations in double precision
    *
    * @param[in] show_progress if true it shows a nice progress bar
    */
    template<typename RealType, typename VectorType>
    template<generation_policy Policy>
    void pair_prod_lookup_table<RealType, VectorType>::generate(
        const bool show_progress)
    {
        constexpr bool use_internal_double =
            (Policy == generation_policy::force_internal_double) &&
            !std::is_same<RealType,double>();

        auto t_start =  std::chrono::system_clock::now();

        const int chi_size = m_params.chi_phot_how_many;
        const int frac_size = m_params.frac_how_many;

        const auto all_coords = get_all_coordinates();
        auto all_vals = std::vector<RealType>(all_coords.size());

        int count = 0;
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < chi_size; ++i){
            const auto chi_phot = all_coords[i*frac_size][0];
            auto chi_parts = std::vector<RealType>(frac_size);
            std::transform(
                all_coords.begin()+i*frac_size,
                all_coords.begin()+(i+1)*frac_size,
                chi_parts.begin(),
                [](auto el){return el[1];}
            );

            std::vector<RealType> vals = std::vector<RealType>(frac_size);
            PXRMP_INTERNAL_CONSTEXPR_IF (use_internal_double){
                vals = aux_generate_double(
                    chi_phot, chi_parts);
            } else {
                vals = compute_cumulative_prob_opt(
                   chi_phot, chi_parts);
            }
            //make sure that the last point is exactly 0.5
            vals.back() = math::half<RealType>;

            std::copy(vals.begin(), vals.end(), all_vals.begin()+i*frac_size);

            #pragma omp critical
            {
                count++;
                utils::draw_progress(count, chi_size, "BW pair prod", 1);
            }
        }

        for (auto& val : all_vals){
            if(std::isnan(val))
                throw std::runtime_error("Error: nan detected in generated table!");
        }

        set_all_vals(all_vals);
        auto t_end =  std::chrono::system_clock::now();
        utils::draw_progress(
            count, chi_size, "BW pair prod", 1, true);
        std::cout << " Done in " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                t_end - t_start).count()/1000.0 << " seconds. \n" << std::endl;

        m_init_flag = true;
    }

    //__________________________________________________________________________
}
}
}
}

#endif //PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE_TABLES_GENERATOR
