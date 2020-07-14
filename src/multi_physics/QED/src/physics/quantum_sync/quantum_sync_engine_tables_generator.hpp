#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR

//This .hpp file extends the implementation of the lookup tables
//for Quantum Synchrotron photon production with methods to generate
//the lookup tables.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "../../qed_commons.h"

//Implements methods of BW lookup tables
#include "quantum_sync_engine_tables.hpp"
//Uses BW tabulated functions
#include "quantum_sync_engine_tabulated_functions.hpp"
//Uses progress bar
#include "../../math/cmath_overloads.hpp"
//Uses progress bar
#include "../../utils/progress_bar.hpp"

#include <omp.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <type_traits>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    //________________ dN/dt table _____________________________________________

    /**
    * Auxiliary function used to compute the G function in double precision
    * with a single precision argument and to cast back the result to single precision
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used (relevant for the class of which is a method is member)
    *
    * @param[in] x the value to be passed to compute_G_function
    *
    * @return the result of compute_G_function
    */
    template<typename RealType, typename VectorType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    RealType dndt_lookup_table<RealType, VectorType>::
    aux_generate_double(RealType x)
    {
        return static_cast<RealType>(
            compute_G_function<double>(x));
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
        for (int i = 0; i < all_coords.size(); ++i){
            PXRMP_INTERNAL_CONSTEXPR_IF (use_internal_double){
                all_vals[i] = aux_generate_double(all_coords[i]);
            }
            else {
                all_vals[i] = compute_G_function(all_coords[i]);
            }

            #pragma omp critical
            {
                count++;
                utils::draw_progress(count,
                    all_vals.size(), "Quantum sync dN/dt", 1);
            }
        }
        set_all_vals(all_vals);

        auto t_end =  std::chrono::system_clock::now();
        utils::draw_progress(
            count, all_vals.size(), "Quantum sync dN/dt", 1, true);
        std::cout << " Done in " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                t_end - t_start).count()/1000.0 << " seconds. \n" << std::endl;

        m_init_flag = true;
    }

    //__________________________________________________________________________

    //________________ Photon emission table ___________________________________

    /**
    * Auxiliary function used to compute cumulative probability distribution in double precision
    * with a single precision argument and to cast back the result to single precision
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used (relevant for the class of which is a method is member)
    *
    * @param[in] x the value to be passed to compute_cumulative_prob
    *
    * @return the result of compute_cumulative_prob
    */
    template<typename RealType, typename VectorType>
    PXRMP_INTERNAL_FORCE_INLINE_DECORATOR
    std::vector<RealType>
    photon_emission_lookup_table<RealType, VectorType>::
    aux_generate_double(const RealType x, const std::vector<RealType>& y)
    {
        auto dtemp = std::vector<double>(y.size());
        std::transform(y.begin(), y.end(), dtemp.begin(),
            [](RealType yy){ return static_cast<double>(yy); });
        const auto dres = compute_cumulative_prob<
            double, std::vector<double>>(x,dtemp);
        auto res = std::vector<RealType>(y.size());
        std::transform(dres.begin(), dres.end(), res.begin(),
            [](double rr){ return static_cast<RealType>(rr); });
        return res;
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
    void photon_emission_lookup_table<RealType, VectorType>::generate(
        const bool show_progress)
    {
        constexpr bool use_internal_double =
            (Policy == generation_policy::force_internal_double) &&
            !std::is_same<RealType,double>();

        auto t_start =  std::chrono::system_clock::now();

        const int chi_size = m_params.chi_part_how_many;
        const int frac_size = m_params.frac_how_many;

        const auto all_coords = get_all_coordinates();
        auto all_vals = std::vector<RealType>(all_coords.size());

        auto fracs = std::vector<RealType>(frac_size);
        for(int j = 0; j < frac_size; ++j){
            fracs[j] = all_coords[j][1];
        }

        int count = 0;
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < chi_size; ++i){
            std::vector<RealType> temp;
            PXRMP_INTERNAL_CONSTEXPR_IF (use_internal_double){
                temp = aux_generate_double(
                    all_coords[i*frac_size][0],fracs);
            } else {
                temp = compute_cumulative_prob(
                    all_coords[i*frac_size][0],fracs);
            }

            std::copy(temp.begin(), temp.end(), all_vals.begin()+i*frac_size);

            #pragma omp critical
            {
                count++;
                utils::draw_progress(count, chi_size, "QS photon emission", 1);
            }
        }

        set_all_vals(all_vals);
        auto t_end =  std::chrono::system_clock::now();
        utils::draw_progress(
            count, chi_size, "QS photon emission", 1, true);
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

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR
