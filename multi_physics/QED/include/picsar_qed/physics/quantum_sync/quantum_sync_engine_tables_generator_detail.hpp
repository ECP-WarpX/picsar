#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR_DETAIL
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR_DETAIL

//This .hpp file contains the actual implementation of lookup tables generation
//for Quantum Synchrotron photon production.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Uses BW tabulated functions
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp"
//Uses cmath overloads
#include "picsar_qed/math/cmath_overloads.hpp"
//Uses progress bar
#include "picsar_qed/utils/progress_bar.hpp"

#ifdef PXRMP_HAS_OPENMP
    #include <omp.h>
#endif

#include <vector>
#include <chrono>
#include <iostream>
#include <type_traits>
#include <stdexcept>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{
namespace detail{
    //________________ dN/dt table _____________________________________________

    /**
    * Implements the generation of the dN/dt lookup table (not usable on GPUs).
    *
    * @tparam RealType the floating point type to be used
    * @tparam ForceInternalDouble if true it forces internal calculations in double precision
    * @tparam ShowProgress if true it shows a nice progress bar
    *
    * @param[in] all_coords a constant reference to the coordinates of the lookup table
    *
    * @return the values corresponding to all_coords
    */
    template<typename RealType, bool ForceInternalDouble, bool ShowProgress>
    std::vector<RealType>
    generate_dndt_lookup_table(
        const std::vector<RealType>& all_coords)
    {
        auto t_start =  std::chrono::system_clock::now();

        auto all_vals = std::vector<RealType>(all_coords.size());

        int count = 0;

#ifdef PXRMP_HAS_OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < static_cast<int>(all_coords.size()); ++i){
            PXRMP_CONSTEXPR_IF (ForceInternalDouble && !std::is_same<RealType,double>()){
                all_vals[i] = static_cast<RealType>(
                    compute_G_function<double>(all_coords[i]));
            }
            else {
                all_vals[i] = compute_G_function(all_coords[i]);
            }

            PXRMP_CONSTEXPR_IF (ShowProgress){
                #pragma omp critical
                {
                    count++;
                    utils::draw_progress(count,
                        all_vals.size(), "Quantum sync dN/dt", 1);
                }
            }
        }

        for (const auto& val : all_vals){
            if(std::isnan(val))
                throw std::runtime_error("Error: nan detected in generated table!");
        }

        auto t_end =  std::chrono::system_clock::now();
        if(ShowProgress){
            utils::draw_progress(
                count, all_vals.size(), "Quantum sync dN/dt", 1, true);

            std::cout << " Done in " <<
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    t_end - t_start).count()/1000.0 << " seconds. \n" << std::endl;
        }

        return all_vals;
    }

    //__________________________________________________________________________

    //________________ Photon emission table ___________________________________

    /**
    * Implements the generation of the photon emission lookup table (not usable on GPUs).
    *
    * @tparam RealType the floating point type to be used
    * @tparam ForceInternalDouble if true it forces internal calculations in double precision
    * @tparam ShowProgress if true it shows a nice progress bar
    *
    * @param[in] all_chi_part a constant reference to the chi coordinates (first axis)
    * @param[in] all_frac a constant reference to the fraction coordinates (second axis)
    *
    * @return the values corresponding to all_chi_part x all_frac
    */
    template<typename RealType, bool ForceInternalDouble, bool ShowProgress>
    std::vector<RealType>
    generate_photon_emission_lookup_table_chipartxfrac(
        const std::vector<RealType>& all_chi_part,
        const std::vector<RealType>& all_frac)
    {
        auto t_start =  std::chrono::system_clock::now();

        const auto all_chi_part_size = all_chi_part.size();
        const auto all_frac_size = all_frac.size();

        auto all_vals = std::vector<RealType>(all_chi_part_size*all_frac_size);

        int count = 0;
#ifdef PXRMP_HAS_OPENMP
        #pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int i = 0; i < static_cast<int>(all_chi_part_size); ++i){
            const auto chi_part = all_chi_part[i];

            auto vals = std::vector<RealType>(all_frac_size);

            PXRMP_CONSTEXPR_IF (ForceInternalDouble && !std::is_same<RealType,double>()){
                const auto d_chi_part = static_cast<double>(chi_part);
                auto d_chi_phots = std::vector<double>(all_frac_size);
                std::transform(
                    all_frac.begin(), all_frac.end(), d_chi_phots.begin(),
                    [=](auto ff){ return static_cast<double>(ff)*d_chi_part;});
                const auto dvals = compute_cumulative_prob_opt<
                    double, std::vector<double>>(d_chi_part, d_chi_phots);
                std::transform(dvals.begin(), dvals.end(), vals.begin(),
                    [](double r){ return static_cast<RealType>(r);});
            } else {
                auto chi_phots = std::vector<RealType>(all_frac_size);
                std::transform(
                    all_frac.begin(), all_frac.end(), chi_phots.begin(),
                    [=](RealType ff){return ff*chi_part;});
                vals = compute_cumulative_prob_opt(
                    chi_part, chi_phots);
            }

            //make sure that the last point is exactly 1.0
            vals.back() = math::one<RealType>;

            std::copy(vals.begin(), vals.end(), all_vals.begin() + i*all_frac_size);

            PXRMP_CONSTEXPR_IF (ShowProgress){
#ifdef PXRMP_HAS_OPENMP
                #pragma omp critical
#endif
                {
                    count++;
                    utils::draw_progress(count, all_chi_part_size, "QS photon emission", 1);
                }
            }
        }

        for (const auto& val : all_vals){
            if(std::isnan(val))
                throw std::runtime_error("Error: nan detected in generated table!");
        }

        auto t_end =  std::chrono::system_clock::now();
        if(ShowProgress){
            utils::draw_progress(
                count, all_chi_part_size, "QS photon emission", 1, true);

            std::cout << " Done in " <<
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    t_end - t_start).count()/1000.0 << " seconds. \n" << std::endl;
        }

        return all_vals;
    }

    //__________________________________________________________________________

}
}
}
}
}

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR
