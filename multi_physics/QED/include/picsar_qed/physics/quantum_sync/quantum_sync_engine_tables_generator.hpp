#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR

//This .hpp file extends the implementation of the lookup tables
//for Quantum Synchrotron photon production with methods to generate
//the lookup tables. In order to avoid code duplication, the actual
//implementation of the generation routines is in
//quantum_sync_engine_tables_generator_detail.hpp .
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Implements methods of BW lookup tables
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp"
//Uses source file containing the actual implementation of table generation methods
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator_detail.hpp"

#ifdef PXRMP_HAS_OPENMP
    #include <omp.h>
#endif

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    //________________ dN/dt table _____________________________________________

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
        const auto& all_coords = get_all_coordinates();

        const auto all_vals = show_progress ?
            detail::generate_dndt_lookup_table<
            RealType, Policy ==  generation_policy::force_internal_double, true>(
                all_coords):
            detail::generate_dndt_lookup_table<
            RealType, Policy ==  generation_policy::force_internal_double, false>(
                all_coords);

        set_all_vals(all_vals);
        m_init_flag = true;
    }

    //__________________________________________________________________________

    //________________ Photon emission table ___________________________________

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
        const auto all_coords = get_all_coordinates();

        auto all_chi_part = std::vector<RealType>(m_params.chi_part_how_many);
        auto all_frac = std::vector<RealType>(m_params.frac_how_many);

        for(int i = 0; i < m_params.chi_part_how_many; ++i)
            all_chi_part[i] = all_coords[i*m_params.frac_how_many][0];

        for(int i = 0; i < m_params.frac_how_many; ++i)
            all_frac[i] = all_coords[i][1];

        const auto all_vals = show_progress ?
            detail::generate_photon_emission_lookup_table_chipartxfrac<
            RealType, Policy ==  generation_policy::force_internal_double, true>(
                all_chi_part, all_frac):
            detail::generate_photon_emission_lookup_table_chipartxfrac<
            RealType, Policy ==  generation_policy::force_internal_double, false>(
                all_chi_part, all_frac);

        set_all_vals(all_vals);
        m_init_flag = true;
    }

    //__________________________________________________________________________


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
    void tailopt_photon_emission_lookup_table<RealType, VectorType>::generate(
        const bool show_progress)
    {
        const auto all_coords = get_all_coordinates();

        auto all_chi_part = std::vector<RealType>(m_params.chi_part_how_many);
        auto all_frac = std::vector<RealType>(m_params.frac_how_many);

        for(int i = 0; i < m_params.chi_part_how_many; ++i)
            all_chi_part[i] = all_coords[i*m_params.frac_how_many][0];

        for(int i = 0; i < m_params.frac_how_many; ++i)
            all_frac[i] = all_coords[i][1];

        const auto all_vals = show_progress ?
            detail::generate_photon_emission_lookup_table_chipartxfrac<
            RealType, Policy ==  generation_policy::force_internal_double, true>(
                all_chi_part, all_frac):
            detail::generate_photon_emission_lookup_table_chipartxfrac<
            RealType, Policy ==  generation_policy::force_internal_double, false>(
                all_chi_part, all_frac);

        set_all_vals(all_vals);
        m_init_flag = true;
    }

    //__________________________________________________________________________

}
}
}
}

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES_GENERATOR
