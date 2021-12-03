#ifndef PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
#define PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES

//This .hpp file contains the implementation of the lookup tables
//for Quantum Synchrotron photon production.
//Please have a look at the jupyter notebook "validation.ipynb"
//in QED_tests/validation for a more in-depth discussion.
//
// References:
// 1) C.P.Ridgers et al. Journal of Computational Physics 260, 1 (2014)
// 2) A.Gonoskov et al. Phys. Rev. E 92, 023305 (2015)

//Should be included by all the src files of the library
#include "picsar_qed/qed_commons.h"

//Uses picsar tables
#include "picsar_qed/containers/picsar_tables.hpp"
//Uses GPU-friendly arrays
#include "picsar_qed/containers/picsar_array.hpp"
//Uses picsar_span
#include "picsar_qed/containers/picsar_span.hpp"
//Uses mathematical constants
#include "picsar_qed/math/math_constants.h"
//Uses serialization
#include "picsar_qed/utils/serialization.hpp"
//Uses interpolation and upper_bound
#include "picsar_qed/utils/picsar_algo.hpp"
//Uses log and exp
#include "picsar_qed/math/cmath_overloads.hpp"
//Uses auxiliary functors
#include "quantum_sync_engine_tables_detail.hpp"

#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace picsar{
namespace multi_physics{
namespace phys{
namespace quantum_sync{

    //________________ Default parameters ______________________________________

    //Reasonable default values for the dndt_lookup_table_params
    //and the pair_prod_lookup_table_params (see below)
    template <typename T>
    constexpr T default_chi_part_min = static_cast<T>(1.0e-3); /*Default minimum particle chi parameter*/
    template <typename T>
    constexpr T default_chi_part_max = static_cast<T>(1.0e3); /* Default maximum particle chi parameter*/
    const int default_chi_part_how_many = 256; /* Default number of grid points for particle chi */
    const int default_frac_how_many = 256; /* Default number of grid points for photon chi */
    template <typename T>
    constexpr T default_frac_min = static_cast<T>(1e-12); /* Default value of the minimum chi_photon fraction */

    template <typename T>
    constexpr T default_frac_switch = static_cast<T>(5e-2); /* Default value of the minimum chi_photon fraction */
    const int default_frac_first = 128; /* Default number of grid points for photon chi */


    //__________________________________________________________________________

    //________________ dN/dt table _____________________________________________


    /**
    * This structure holds the parameters to generate a dN/dt
    * lookup table. The dN/dt lookup table stores the values of the
    * T function (see validation script)
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    struct dndt_lookup_table_params{
        RealType chi_part_min = static_cast<RealType>(0.0); /*Minimum particle chi parameter*/
        RealType chi_part_max = static_cast<RealType>(0.0); /*Maximum particle chi parameter*/
        int chi_part_how_many = 0; /* Number of grid points for particle chi */

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

    /**
    * generation_policy::force_internal_double can be used to force the
    * calculations of a lookup tables using double precision, even if
    * the final result is stored in a single precision variable.
    */
    enum class generation_policy{
        regular,
        force_internal_double
    };

    /**
    * This class provides the lookup table for dN/dt,
    * storing the values of the G function (see validation script)
    * and providing methods to perform interpolations.
    * It also provides methods for serialization (export to byte array,
    * import from byte array) and to generate "table views" based on
    * non-owning pointers (this is crucial in order to use the table
    * in GPU kernels, as explained below).
    *
    * Internally, this table stores log(G(log(chi))).
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
            * The use of views is crucial for GPUs. As demonstrated
            * in test_gpu/test_quantum_sync.cu, it is possible to:
            * - define a thin wrapper around a thrust::device_vector
            * - use this wrapper as a template parameter for a dndt_lookup_table
            * - initialize the lookup table on the CPU
            * - generate a view
            * - pass by copy this view to a GPU kernel (see get_view() method)
            *
            * @tparam RealType the floating point type to be used
            */
            typedef dndt_lookup_table<
                RealType, containers::picsar_span<const RealType>> view_type;

            /**
            * Empty constructor
            **/
            constexpr
            dndt_lookup_table(){}

            /**
            * Constructor (not designed for GPU usage)
            * After construction the table is uninitialized. The user has to generate
            * the G function values before being able to use the table.
            *
            * @param params table parameters
            */
            dndt_lookup_table(dndt_lookup_table_params<RealType> params):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    math::m_log(params.chi_part_min),
                    math::m_log(params.chi_part_max),
                    VectorType(params.chi_part_how_many)}}
            {}

            /**
            * Constructor (not designed for GPU usage)
            * This constructor allows the user to initialize the table with
            * a vector of values.
            *
            * @param params parameters for table generation
            * @param vals values of the G function
            */
            dndt_lookup_table(dndt_lookup_table_params<RealType> params,
                VectorType vals):
            m_params{params},
            m_table{containers::equispaced_1d_table<RealType, VectorType>{
                    math::m_log(params.chi_part_min),
                    math::m_log(params.chi_part_max),
                    vals}}
            {
                m_init_flag = true;
            }

            /*
            * Generates the content of the lookup table (not usable on GPUs).
            * This function is implemented elsewhere
            * (in breit_wheeler_engine_tables_generator.hpp)
            * since it requires a recent version of the Boost library.
            *
            * @tparam Policy the generation policy (can force calculations in double precision)
            *
            * @param[in] show_progress if true a progress bar is shown
            */
            template <generation_policy Policy = generation_policy::regular>
            void generate(const bool show_progress  = true);

            /*
            * Initializes the lookup table from a byte array.
            * This method is not usable on GPUs.
            *
            * @param[in] raw_data the byte array
            */
            dndt_lookup_table(const std::vector<char>& raw_data)
            {
                using namespace utils;

                constexpr size_t min_size =
                    sizeof(char)+ //single or double precision
                    sizeof(m_params);

                if (raw_data.size() < min_size)
                    throw std::runtime_error("Binary data is too small \
                    to be a Quantum Synchrotron G-function lookup-table.");

                auto it_raw_data = raw_data.begin();

                if (serialization::get_out<char>(it_raw_data) !=
                    static_cast<char>(sizeof(RealType))){
                    throw std::runtime_error("Mismatch between RealType used \
                    to write and to read the Quantum Synchrotron G-function lookup-table");
                }

                m_params = serialization::get_out<
                    dndt_lookup_table_params<RealType>>(it_raw_data);
                m_table = containers::equispaced_1d_table<
                    RealType, VectorType>{std::vector<char>(it_raw_data,
                        raw_data.end())};

                m_init_flag = true;
            }

            /**
            * Operator==
            *
            * @param[in] rhs a structure of the same type
            *
            * @return true if rhs is equal to *this. false otherwise
            */
            PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
            bool operator== (
                const dndt_lookup_table<RealType, VectorType> &rhs) const
            {
                return
                    (m_params == rhs.m_params) &&
                    (m_init_flag == rhs.m_init_flag) &&
                    (m_table == rhs.m_table);
            }

            /*
            * Returns a table view for the current table
            * (i.e. a table built using non-owning picsar_span
            * vectors). A view_type is very lightweight and can
            * be passed by copy to functions and GPU kernels.
            * Indeed it contains non-owning pointers to the data
            * held by the original table.
            * The method is not designed to be run on GPUs.
            *
            * @return a table view
            */
            view_type get_view() const
            {
                if(!m_init_flag)
                    throw std::runtime_error("Can't generate a view of an \
                    uninitialized table");
                const auto span = containers::picsar_span<const RealType>{
                    static_cast<size_t>(m_params.chi_part_how_many),
                    m_table.get_values_reference().data()
                };
                const view_type view{m_params, span};
                return view;
            }

            /*
            * Uses the lookup table to interpolate G function
            * at a given position chi_phot. If chi_part is out
            * of table either the minimum or the maximum value
            * is used. In addition, it checks if chi_phot is out of table
            * and stores the result in a bool variable.
            *
            * @param[in] chi_part where the G function is interpolated
            * @param[out] is_out set to true if chi_part is out of table
            *
            * @return the value of the G function
            */
            PXRMP_GPU_QUALIFIER
            PXRMP_FORCE_INLINE
            RealType interp(
                RealType chi_part, bool* const is_out = nullptr) const noexcept
            {
                if(chi_part<m_params.chi_part_min){
                    chi_part = m_params.chi_part_min;
                    if (is_out != nullptr) *is_out = true;
                }
                else if (chi_part > m_params.chi_part_max){
                    chi_part = m_params.chi_part_max;
                    if (is_out != nullptr) *is_out = true;
                }
                return math::m_exp(m_table.interp(math::m_log(chi_part)));
            }

            /**
            * Exports all the coordinates (chi_particle) of the table to a std::vector
            * (not usable on GPUs)
            *
            * @return a vector containing all the table coordinates
            */
            std::vector<RealType> get_all_coordinates() const noexcept
            {
                auto all_coords = m_table.get_all_coordinates();
                std::transform(all_coords.begin(),all_coords.end(),all_coords.begin(),
                    [](RealType a){return math::m_exp(a);});
                return all_coords;
            }

            /**
            * Imports table values from an std::vector. Values
            * should correspond to coordinates exported with
            * get_all_coordinates(). Not usable on GPU.
            *
            * @param[in] a std::vector containing table values
            *
            * @return false if the value vector has the wrong length. True otherwise.
            */
            bool set_all_vals(const std::vector<RealType>& vals)
            {
                if(static_cast<int>(vals.size()) == m_table.get_how_many_x()){
                    for(int i = 0; i < static_cast<int>(vals.size()); ++i){
                        m_table.set_val(i, math::m_log(vals[i]));
                    }
                    m_init_flag = true;
                    return true;
                }
                return false;
            }

            /*
            * Checks if the table has been initialized.
            *
            * @return true if the table has been initialized, false otherwise
            */
            PXRMP_GPU_QUALIFIER
            PXRMP_FORCE_INLINE
            bool is_init() const
            {
                return m_init_flag;
            }

            /*
            * Converts the table to a byte vector
            *
            * @return a byte vector
            */
            std::vector<char> serialize() const
            {
                using namespace utils;

                if(!m_init_flag)
                    throw std::runtime_error("Cannot serialize \
                    an uninitialized table");

                std::vector<char> res;

                serialization::put_in(static_cast<char>(sizeof(RealType)), res);
                serialization::put_in(m_params, res);

                auto tdata = m_table.serialize();
                res.insert(res.end(), tdata.begin(), tdata.end());

                return res;
            }

        protected:
            dndt_lookup_table_params<RealType> m_params; /* Table parameters*/
            bool m_init_flag = false;  /* Initialization flag*/
            containers::equispaced_1d_table<
                RealType, VectorType> m_table; /* Table data */

        private:
            /*
            * Auxiliary function used for the generation of the lookup table.
            * (not usable on GPUs). This function is implemented elsewhere
            * (in quantum_sync_engine_tables_generator.hpp)
            * since it requires the Boost library.
            */
            PXRMP_FORCE_INLINE
            static RealType aux_generate_double(RealType x);

    };

    //__________________________________________________________________________

    //________________ Photon emission table ___________________________________

    /**
    * This structure holds the parameters to generate a photon
    * emission lookup table. The lookup table stores the
    * values of a cumulative probability distribution (see validation script)
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    struct photon_emission_lookup_table_params{
        RealType chi_part_min = static_cast<RealType>(0.0); /*Minimum particle chi parameter*/
        RealType chi_part_max = static_cast<RealType>(0.0); /*Maximum particle chi parameter*/
        RealType frac_min = static_cast<RealType>(0.0); /*Minimum chi photon fraction to be stored*/
        int chi_part_how_many = 0; /* Number of grid points for particle chi */
        int frac_how_many = 0; /* Number of grid points for photon chi fraction */

        /**
        * Operator==
        *
        * @param[in] rhs a structure of the same type
        * @return true if rhs is equal to *this. false otherwise
        */
        bool operator== (
            const photon_emission_lookup_table_params<RealType> &rhs) const
        {
            return (chi_part_min == rhs.chi_part_min) &&
                (chi_part_max == rhs.chi_part_max) &&
                (frac_min == rhs.frac_min) &&
                (chi_part_how_many == rhs.chi_part_how_many) &&
                (frac_how_many == rhs.frac_how_many);
        }
    };

    /**
    * The default photon_emission_lookup_table_params
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    constexpr auto default_photon_emission_lookup_table_params =
        photon_emission_lookup_table_params<RealType>{default_chi_part_min<RealType>,
                                           default_chi_part_max<RealType>,
                                           default_frac_min<RealType>,
                                           default_chi_part_how_many,
                                           default_frac_how_many};

    /**
    * This class provides the lookup table for photon emission
    * storing the values of a cumulative probability distribution
    * (see validation script) and providing methods to perform interpolations.
    * It also provides methods for serialization (export to byte array,
    * import from byte array) and to generate "table views" based on
    * non-owning pointers (this is crucial in order to use the table
    * in GPU kernels, as explained below).
    *
    * Internally, this table stores
    * log(P(log(chi_particle), log(chi_photon/chi_particle))).
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used internally (e.g. std::vector)
    */
    template<typename RealType, typename VectorType>
    class photon_emission_lookup_table
    {

        public:

            /**
            * A view_type is essentially a photon_emission_lookup_table which
            * uses non-owning, constant, pointers to hold the data.
            * The use of views is crucial for GPUs. As demonstrated
            * in test_gpu/test_quantum_sync.cu, it is possible to:
            * - define a thin wrapper around a thrust::device_vector
            * - use this wrapper as a template parameter for a photon_emission_lookup_table
            * - initialize the lookup table on the CPU
            * - generate a view
            * - pass by copy this view to a GPU kernel (see get_view() method)
            *
            * @tparam RealType the floating point type to be used
            */
            typedef photon_emission_lookup_table<
                RealType, containers::picsar_span<const RealType>> view_type;

            /**
            * Empty constructor
            */
            constexpr
            photon_emission_lookup_table(){}

            /**
            * Constructor (not designed for GPU usage)
            * After construction the table is uninitialized. The user has to generate
            * the cumulative probability distribution before being able to use the table.
            *
            * @param params table parameters
            */
            photon_emission_lookup_table(
                photon_emission_lookup_table_params<RealType> params):
                m_params{params},
                m_table{containers::equispaced_2d_table<RealType, VectorType>{
                    math::m_log(params.chi_part_min),
                    math::m_log(params.chi_part_max),
                    math::m_log(params.frac_min),
                    math::m_log(math::one<RealType>),
                    params.chi_part_how_many, params.frac_how_many,
                    VectorType(params.chi_part_how_many * params.frac_how_many)}}
                {}

            /**
            * Constructor (not designed for GPU usage)
            * This constructor allows the user to initialize the table with
            * a vector of values.
            *
            * @param params table parameters
            */
            photon_emission_lookup_table(
                photon_emission_lookup_table_params<RealType> params,
                VectorType vals):
                m_params{params},
                m_table{containers::equispaced_2d_table<RealType, VectorType>{
                    math::m_log(params.chi_part_min),
                    math::m_log(params.chi_part_max),
                    math::m_log(params.frac_min),
                    math::m_log(math::one<RealType>),
                    params.chi_part_how_many, params.frac_how_many,
                    vals}}
            {
                m_init_flag = true;
            }

            /*
            * Generates the content of the lookup table (not usable on GPUs).
            * This function is implemented elsewhere
            * (in breit_wheeler_engine_tables_generator.hpp)
            * since it requires a recent version of the Boost library.
            *
            * @tparam Policy the generation policy (can force calculations in double precision)
            *
            * @param[in] show_progress if true a progress bar is shown
            */
            template <generation_policy Policy = generation_policy::regular>
            void generate(bool show_progress = true);

            /*
            * Initializes the lookup table from a byte array.
            * This method is not usable on GPUs.
            *
            * @param[in] raw_data the byte array
            */
            photon_emission_lookup_table(const std::vector<char>& raw_data)
            {
                using namespace utils;

                constexpr size_t min_size =
                    sizeof(char)+//single or double precision
                    sizeof(m_params);

                if (raw_data.size() < min_size)
                    throw std::runtime_error("Binary data is too small \
                    to be a Quantum Synchrotron emisson lookup-table.");

                auto it_raw_data = raw_data.begin();

                if (serialization::get_out<char>(it_raw_data) !=
                    static_cast<char>(sizeof(RealType))){
                    throw std::runtime_error("Mismatch between RealType \
                    used to write and to read the Quantum Synchrotron lookup-table");
                }

                m_params = serialization::get_out<
                    photon_emission_lookup_table_params<RealType>>(it_raw_data);
                m_table = containers::equispaced_2d_table<
                    RealType, VectorType>{std::vector<char>(it_raw_data,
                    raw_data.end())};

                m_init_flag = true;
            }

            /**
            * Operator==
            *
            * @param[in] rhs a structure of the same type
            *
            * @return true if rhs is equal to *this. false otherwise
            */
            PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
            bool operator== (
                const photon_emission_lookup_table<
                    RealType, VectorType> &rhs) const
            {
                return
                    (m_params == rhs.m_params) &&
                    (m_init_flag == rhs.m_init_flag) &&
                    (m_table == rhs.m_table);
            }

            /*
            * Returns a table view for the current table
            * (i.e. a table built using non-owning picsar_span
            * vectors). A view_type is very lightweight and can
            * be passed by copy to functions and GPU kernels.
            * Indeed it contains non-owning pointers to the data
            * held by the original table.
            * The method is not designed to be run on GPUs.
            *
            * @return a table view
            */
            view_type get_view() const
            {
                if(!m_init_flag)
                    throw std::runtime_error("Can't generate a view of an \
                    uninitialized table");
                const auto span = containers::picsar_span<const RealType>{
                    static_cast<size_t>(m_params.chi_part_how_many *
                        m_params.frac_how_many),
                        m_table.get_values_reference().data()};

                return view_type{m_params, span};
            }

            /*
            * Uses the lookup table to extract the chi value of
            * the generated photon from a cumulative probability
            * distribution, given the chi parameter of the particle and a
            * random number uniformly distributed in [0,1). If chi_part is out
            * of table either the minimum or the maximum value is used.
            * The method uses the lookup table to invert the equation:
            * unf_zero_one_minus_epsi = P(chi_part, X)
            * where X is the ratio between chi_photon and chi_particle.
            * If X is out of table, 0 is returned (i.e. a photon with
            * very low energy is emitted, so it can be disregarded)
            * The method also checks if chi_phot is out of table
            * and stores the result in a bool variable.
            *
            * @param[in] chi_phot to be used for interpolation
            * @param[in] unf_zero_one_minus_epsi a uniformly distributed random number in [0,1)
            * @param[out] is_out set to true if chi_part is out of table
            *
            * @return chi of one of the generated particles
            */
            PXRMP_GPU_QUALIFIER
            PXRMP_FORCE_INLINE
            RealType interp(
                const RealType chi_part,
                const RealType unf_zero_one_minus_epsi,
                bool* const is_out = nullptr) const noexcept
            {
                using namespace math;

                auto e_chi_part = chi_part;
                if(chi_part<m_params.chi_part_min){
                    e_chi_part = m_params.chi_part_min;
                    if (is_out != nullptr) *is_out = true;
                }
                else if (chi_part > m_params.chi_part_max){
                    e_chi_part = m_params.chi_part_max;
                    if (is_out != nullptr) *is_out = true;
                }

                const auto log_e_chi_part = m_log(e_chi_part);
                const auto log_prob = m_log(one<RealType>-unf_zero_one_minus_epsi);

                const auto upper_frac_index = utils::picsar_lower_bound_functor(
                    0, m_params.frac_how_many,log_prob,[&](int i){
                        return (m_table.interp_first_coord(
                            log_e_chi_part, i));
                        });

                if(upper_frac_index == 0)
                    return zero<RealType>;

                if(upper_frac_index ==  m_params.frac_how_many)
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

            /**
            * Exports all the coordinates (chi_particle, chi_photon)
            * of the table to a std::vector
            * of 2-elements arrays (not usable on GPUs).
            *
            * @return a vector containing all the table coordinates
            */
            std::vector<std::array<RealType,2>> get_all_coordinates() const noexcept
            {
                auto all_coords = m_table.get_all_coordinates();
                std::transform(all_coords.begin(),all_coords.end(),all_coords.begin(),
                    [](std::array<RealType,2> a){return
                        std::array<RealType,2>{math::m_exp(a[0]), math::m_exp(a[1]) * math::m_exp(a[0]) };});
                return all_coords;
            }

            /**
            * Imports table values from an std::vector. Values
            * should correspond to coordinates exported with
            * get_all_coordinates(). Not usable on GPU.
            *
            * @param[in] a std::vector containing table values
            *
            * @return false if the value vector has the wrong length. True otherwise.
            */
            bool set_all_vals(const std::vector<RealType>& vals)
            {
                if(static_cast<int>(vals.size()) == m_table.get_how_many_x()*
                    m_table.get_how_many_y()){
                    for(int i = 0; i < static_cast<int>(vals.size()); ++i){
                        auto val = math::m_log(vals[i]);
                        if(std::isinf(val))
                            val = std::numeric_limits<RealType>::lowest();
                        m_table.set_val(i, val);
                    }
                    m_init_flag = true;
                    return true;
                }
                return false;
            }

            /*
            * Checks if the table has been initialized.
            *
            * @return true if the table has been initialized, false otherwise
            */
            PXRMP_GPU_QUALIFIER
            PXRMP_FORCE_INLINE
            bool is_init() const
            {
                return m_init_flag;
            }

            /*
            * Converts the table to a byte vector
            *
            * @return a byte vector
            */
            std::vector<char> serialize() const
            {
                using namespace utils;

                if(!m_init_flag)
                    throw std::runtime_error("Cannot serialize an unitialized table");

                std::vector<char> res;

                serialization::put_in(static_cast<char>(sizeof(RealType)), res);
                serialization::put_in(m_params, res);

                auto tdata = m_table.serialize();
                res.insert(res.end(), tdata.begin(), tdata.end());

                return res;
            }

        protected:
            photon_emission_lookup_table_params<RealType> m_params; /* Table parameters*/
            bool m_init_flag = false; /* Initialization flag*/
            containers::equispaced_2d_table<
                RealType, VectorType> m_table; /* Table data*/

        private:
            /*
            * Auxiliary function used for the generation of the lookup table.
            * (not usable on GPUs). This function is implemented elsewhere
            * (in breit_wheeler_engine_tables_generator.hpp)
            * since it requires the Boost library.
            */
            PXRMP_FORCE_INLINE
            static std::vector<RealType>
            aux_generate_double(RealType x,
                const std::vector<RealType>& y);

    };

    //__________________________________________________________________________

    //________________ Tail-optimized photon emission table ____________________

    /**
    * This structure holds the parameters to generate a photon
    * emission lookup table. The lookup table stores the
    * values of a cumulative probability distribution (see validation script)
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    struct tailopt_photon_emission_lookup_table_params{
        RealType chi_part_min = static_cast<RealType>(0.0); /*Minimum particle chi parameter*/
        RealType chi_part_max = static_cast<RealType>(0.0); /*Maximum particle chi parameter*/
        RealType frac_min = static_cast<RealType>(0.0); /*Minimum chi photon fraction to be stored*/
        RealType frac_switch = static_cast<RealType>(0.0);
        int chi_part_how_many = 0; /* Number of grid points for particle chi */
        int frac_how_many = 0; /* Number of grid points for photon chi fraction */
        int frac_first = 0;

        /**
        * Operator==
        *
        * @param[in] rhs a structure of the same type
        * @return true if rhs is equal to *this. false otherwise
        */
        bool operator== (
            const tailopt_photon_emission_lookup_table_params<RealType> &rhs) const
        {
            return (chi_part_min == rhs.chi_part_min) &&
                (chi_part_max == rhs.chi_part_max) &&
                (frac_min == rhs.frac_min) &&
                (frac_switch == rhs.frac_switch) &&
                (chi_part_how_many == rhs.chi_part_how_many) &&
                (frac_how_many == rhs.frac_how_many) &&
                (frac_first == rhs.frac_first);
        }
    };

    /**
    * The default photon_emission_lookup_table_params
    *
    * @tparam RealType the floating point type to be used
    */
    template<typename RealType>
    constexpr auto default_tailopt_photon_emission_lookup_table_params =
        tailopt_photon_emission_lookup_table_params<RealType>{
                                        default_chi_part_min<RealType>,
                                        default_chi_part_max<RealType>,
                                        default_frac_min<RealType>,
                                        default_frac_switch<RealType>,
                                        default_chi_part_how_many,
                                        default_frac_how_many,
                                        default_frac_first};


    template<typename RealType, typename VectorType>
    using Generic2DTableType =
        containers::generic_2d_table<
            RealType, VectorType,
            detail::LinFunctor<RealType>, detail::TailOptFunctor<RealType>,
            detail::ILinFunctor<RealType>, detail::ITailOptFunctor<RealType> >;


    /**
    * This class provides the lookup table for photon emission
    * storing the values of a cumulative probability distribution
    * (see validation script) and providing methods to perform interpolations.
    * It also provides methods for serialization (export to byte array,
    * import from byte array) and to generate "table views" based on
    * non-owning pointers (this is crucial in order to use the table
    * in GPU kernels, as explained below).
    *
    * Internally, this table stores
    * log(P(log(chi_particle), log(chi_photon/chi_particle))).
    *
    * @tparam RealType the floating point type to be used
    * @tparam VectorType the vector type to be used internally (e.g. std::vector)
    */
    template<typename RealType, typename VectorType>
    class tailopt_photon_emission_lookup_table
    {
        public:

            /**
            * A view_type is essentially a photon_emission_lookup_table which
            * uses non-owning, constant, pointers to hold the data.
            * The use of views is crucial for GPUs. As demonstrated
            * in test_gpu/test_quantum_sync.cu, it is possible to:
            * - define a thin wrapper around a thrust::device_vector
            * - use this wrapper as a template parameter for a photon_emission_lookup_table
            * - initialize the lookup table on the CPU
            * - generate a view
            * - pass by copy this view to a GPU kernel (see get_view() method)
            *
            * @tparam RealType the floating point type to be used
            */
            typedef tailopt_photon_emission_lookup_table<
                RealType, containers::picsar_span<const RealType>> view_type;

            /**
            * Empty constructor
            */
            constexpr
            tailopt_photon_emission_lookup_table(){}

            /**
            * Constructor (not designed for GPU usage)
            * After construction the table is uninitialized. The user has to generate
            * the cumulative probability distribution before being able to use the table.
            *
            * @param params table parameters
            */
            tailopt_photon_emission_lookup_table(
                tailopt_photon_emission_lookup_table_params<RealType> params):
                m_params{params},
                m_table{Generic2DTableType<RealType, VectorType>{
                    params.chi_part_how_many, params.frac_how_many,
                    VectorType(params.chi_part_how_many * params.frac_how_many),
                    detail::LinFunctor<RealType>(
                        m_params.chi_part_how_many,
                        math::m_log(m_params.chi_part_min),
                        math::m_log(m_params.chi_part_max)),
                    detail::TailOptFunctor<RealType>(
                        m_params.frac_how_many, m_params.frac_first,
                        m_params.frac_min, math::one<RealType>, m_params.frac_switch),
                    detail::ILinFunctor<RealType>(
                        m_params.chi_part_how_many,
                        math::m_log(m_params.chi_part_min),
                        math::m_log(m_params.chi_part_max)),
                    detail::ITailOptFunctor<RealType>(
                        m_params.frac_how_many, m_params.frac_first,
                        m_params.frac_min, math::one<RealType>, m_params.frac_switch)}}
                {}

            /**
            * Constructor (not designed for GPU usage)
            * This constructor allows the user to initialize the table with
            * a vector of values.
            *
            * @param params table parameters
            */
            tailopt_photon_emission_lookup_table(
                tailopt_photon_emission_lookup_table_params<RealType> params,
                VectorType vals):
                m_params{params},
                m_table{Generic2DTableType<RealType, VectorType>{
                    m_params.chi_part_how_many, m_params.frac_how_many, vals,
                    detail::LinFunctor<RealType>(
                        m_params.chi_part_how_many,
                        math::m_log(m_params.chi_part_min),
                        math::m_log(m_params.chi_part_max)),
                    detail::TailOptFunctor<RealType>(
                        m_params.frac_how_many, m_params.frac_first,
                        m_params.frac_min, math::one<RealType>, m_params.frac_switch),
                    detail::ILinFunctor<RealType>(
                        m_params.chi_part_how_many,
                        math::m_log(m_params.chi_part_min),
                        math::m_log(m_params.chi_part_max)),
                    detail::ITailOptFunctor<RealType>(
                        m_params.frac_how_many, m_params.frac_first,
                        m_params.frac_min, math::one<RealType>, m_params.frac_switch)}}
            {
                m_init_flag = true;
            }

            /*
            * Generates the content of the lookup table (not usable on GPUs).
            * This function is implemented elsewhere
            * (in breit_wheeler_engine_tables_generator.hpp)
            * since it requires a recent version of the Boost library.
            *
            * @tparam Policy the generation policy (can force calculations in double precision)
            *
            * @param[in] show_progress if true a progress bar is shown
            */
            template <generation_policy Policy = generation_policy::regular>
            void generate(bool show_progress = true);

            /*
            * Initializes the lookup table from a byte array.
            * This method is not usable on GPUs.
            *
            * @param[in] raw_data the byte array
            */
            tailopt_photon_emission_lookup_table(const std::vector<char>& raw_data)
            {
                using namespace utils;

                constexpr size_t min_size =
                    sizeof(char)+//single or double precision
                    sizeof(m_params);

                if (raw_data.size() < min_size)
                    throw std::runtime_error("Binary data is too small \
                    to be a Quantum Synchrotron emisson lookup-table.");

                auto it_raw_data = raw_data.begin();

                if (serialization::get_out<char>(it_raw_data) !=
                    static_cast<char>(sizeof(RealType))){
                    throw std::runtime_error("Mismatch between RealType \
                    used to write and to read the Quantum Synchrotron lookup-table");
                }

                m_params = serialization::get_out<
                    tailopt_photon_emission_lookup_table_params<RealType>>(it_raw_data);
                m_table = Generic2DTableType<RealType, VectorType>{
                    std::vector<char>(it_raw_data, raw_data.end())};

                m_init_flag = true;
            }

            /**
            * Operator==
            *
            * @param[in] rhs a structure of the same type
            *
            * @return true if rhs is equal to *this. false otherwise
            */
            PXRMP_GPU_QUALIFIER PXRMP_FORCE_INLINE
            bool operator== (
                const tailopt_photon_emission_lookup_table<
                    RealType, VectorType> &rhs) const
            {
                return
                    (m_params == rhs.m_params) &&
                    (m_init_flag == rhs.m_init_flag) &&
                    (m_table == rhs.m_table);
            }

            /*
            * Returns a table view for the current table
            * (i.e. a table built using non-owning picsar_span
            * vectors). A view_type is very lightweight and can
            * be passed by copy to functions and GPU kernels.
            * Indeed it contains non-owning pointers to the data
            * held by the original table.
            * The method is not designed to be run on GPUs.
            *
            * @return a table view
            */
            view_type get_view() const
            {
                if(!m_init_flag)
                    throw std::runtime_error("Can't generate a view of an \
                    uninitialized table");
                const auto span = containers::picsar_span<const RealType>{
                    static_cast<size_t>(m_params.chi_part_how_many *
                        m_params.frac_how_many),
                        m_table.get_values_reference().data()};

                return view_type{m_params, span};
            }

            /*
            * Uses the lookup table to extract the chi value of
            * the generated photon from a cumulative probability
            * distribution, given the chi parameter of the particle and a
            * random number uniformly distributed in [0,1). If chi_part is out
            * of table either the minimum or the maximum value is used.
            * The method uses the lookup table to invert the equation:
            * unf_zero_one_minus_epsi = P(chi_part, X)
            * where X is the ratio between chi_photon and chi_particle.
            * If X is out of table, 0 is returned (i.e. a photon with
            * very low energy is emitted, so it can be disregarded)
            * The method also checks if chi_phot is out of table
            * and stores the result in a bool variable.
            *
            * @param[in] chi_phot to be used for interpolation
            * @param[in] unf_zero_one_minus_epsi a uniformly distributed random number in [0,1)
            * @param[out] is_out set to true if chi_part is out of table
            *
            * @return chi of one of the generated particles
            */
            PXRMP_GPU_QUALIFIER
            PXRMP_FORCE_INLINE
            RealType interp(
                const RealType chi_part,
                const RealType unf_zero_one_minus_epsi,
                bool* const is_out = nullptr) const noexcept
            {
                using namespace math;

                auto e_chi_part = chi_part;
                if(chi_part<m_params.chi_part_min){
                    e_chi_part = m_params.chi_part_min;
                    if (is_out != nullptr) *is_out = true;
                }
                else if (chi_part > m_params.chi_part_max){
                    e_chi_part = m_params.chi_part_max;
                    if (is_out != nullptr) *is_out = true;
                }

                const auto log_e_chi_part = m_log(e_chi_part);
                const auto log_prob = m_log(one<RealType>-unf_zero_one_minus_epsi);

                const auto upper_frac_index = utils::picsar_upper_bound_functor(
                    0, m_params.frac_how_many,log_prob,[&](int i){
                        return (m_table.interp_first_coord(
                            log_e_chi_part, i));
                        });

                if(upper_frac_index == 0)
                    return zero<RealType>;

                if(upper_frac_index ==  m_params.frac_how_many)
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

            /**
            * Exports all the coordinates (chi_particle, chi_photon)
            * of the table to a std::vector
            * of 2-elements arrays (not usable on GPUs).
            *
            * @return a vector containing all the table coordinates
            */
            std::vector<std::array<RealType,2>> get_all_coordinates() const noexcept
            {
                auto all_coords = m_table.get_all_coordinates();
                std::transform(all_coords.begin(),all_coords.end(),all_coords.begin(),
                    [](std::array<RealType,2> a){return
                        std::array<RealType,2>{math::m_exp(a[0]), math::m_exp(a[1]) * math::m_exp(a[0]) };});
                return all_coords;
            }

            /**
            * Imports table values from an std::vector. Values
            * should correspond to coordinates exported with
            * get_all_coordinates(). Not usable on GPU.
            *
            * @param[in] a std::vector containing table values
            *
            * @return false if the value vector has the wrong length. True otherwise.
            */
            bool set_all_vals(const std::vector<RealType>& vals)
            {
                if(static_cast<int>(vals.size()) == m_table.get_how_many_x()*
                    m_table.get_how_many_y()){
                    for(int i = 0; i < static_cast<int>(vals.size()); ++i){
                        auto val = math::m_log(vals[i]);
                        if(std::isinf(val))
                            val = std::numeric_limits<RealType>::lowest();
                        m_table.set_val(i, val);
                    }
                    m_init_flag = true;
                    return true;
                }
                return false;
            }

            /*
            * Checks if the table has been initialized.
            *
            * @return true if the table has been initialized, false otherwise
            */
            PXRMP_GPU_QUALIFIER
            PXRMP_FORCE_INLINE
            bool is_init() const
            {
                return m_init_flag;
            }

            /*
            * Converts the table to a byte vector
            *
            * @return a byte vector
            */
            std::vector<char> serialize() const
            {
                using namespace utils;

                if(!m_init_flag)
                    throw std::runtime_error("Cannot serialize an unitialized table");

                std::vector<char> res;

                serialization::put_in(static_cast<char>(sizeof(RealType)), res);
                serialization::put_in(m_params, res);

                auto tdata = m_table.serialize();
                res.insert(res.end(), tdata.begin(), tdata.end());

                return res;
            }

        protected:
            tailopt_photon_emission_lookup_table_params<RealType> m_params; /* Table parameters*/
            bool m_init_flag = false; /* Initialization flag*/
            Generic2DTableType<RealType, VectorType> m_table; /* Table data*/

        private:
            /*
            * Auxiliary function used for the generation of the lookup table.
            * (not usable on GPUs). This function is implemented elsewhere
            * (in breit_wheeler_engine_tables_generator.hpp)
            * since it requires the Boost library.
            */
            PXRMP_FORCE_INLINE
            static std::vector<RealType>
            aux_generate_double(RealType x,
                const std::vector<RealType>& y);

        };

        //______________________________________________________________________



}
}
}
}

#endif // PICSAR_MULTIPHYSICS_QUANTUM_SYNC_ENGINE_TABLES
