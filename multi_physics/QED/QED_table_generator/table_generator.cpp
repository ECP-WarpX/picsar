/**
* This program generates the lookup tables of the QED library and saves them both
* in binary format and as a csv file which can be inspected with pyplot
* or gnuplot.
*/

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#ifdef PXRMP_TABLE_GEN_HAS_OPENMP
    #include <omp.h>
#endif

#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp"
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"
#include "picsar_qed/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp"
#include "picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"

namespace px_bw = picsar::multi_physics::phys::breit_wheeler;
namespace px_qs = picsar::multi_physics::phys::quantum_sync;
namespace px_ut = picsar::multi_physics::utils;

// These string constants are used to parse command line instructions
const std::string CMD_HELP_S = "-h";
const std::string CMD_HELP_L = "--help";
const std::string CMD_TABLE = "--table_type";
const std::string CMD_PRECISION = "--precision";
const std::string CMD_CHI_SIZE = "--chi_size";
const std::string CMD_FRAC_SIZE = "--frac_size";
const std::string CMD_CHI_MIN = "--chi_min";
const std::string CMD_CHI_MAX = "--chi_max";
const std::string CMD_FRAC_MIN = "--frac_min";
const std::string CMD_FILENAME = "--name";
const std::string OPT_TABLE_BREIT_WHEELER = "breit_wheeler";
const std::string OPT_TABLE_QUANTUM_SYNC = "quantum_synchrotron";
const std::string OPT_PRECISION_DOUBLE = "double";
const std::string OPT_PRECISION_SINGLE = "single";
const std::string OPT_PRECISION_SINGLE_WITH_DOUBLE_COMP = "single_with_double_comp";
//___________________________________________________________________________________

enum Precision {double_precision, single_precision, single_prec_out_double_prec_comp};
enum TableType {breit_wheeler_table, quantum_synchrotron_table};

template<typename RealType>
struct BreitWheelerTableParams{
    int chi_size = 0;
    int frac_size = 0;
    RealType chi_min = RealType{0.0};
    RealType chi_max = RealType{0.0};
};

template<typename RealType>
struct QuantumSyncTableParams{
    int chi_size = 0;
    int frac_size = 0;
    RealType chi_min = RealType{0.0};
    RealType chi_max = RealType{0.0};
    RealType frac_min = RealType{0.0};
};

/**
* Prints an error message
*
* @param[in] str the error message
*/
void print_error(const std::string& str)
{
    std::cout << "\n [ERROR!] " << str << "\n" << std::endl;
}


/**
* Writes a 1D lookup table in csv format
*
* @tparam RealType the floating point type to be used
* @tparam TableType the table type to be used
* @param[in] left the left limit to sample the table
* @param[in] right the right limit to sample the table
* @param[in] how_many controls how many points are used to sample the lookup table
* @param[in] log_scale controls if the axis of the table is in log or linear scale
* @param[in] file_name the name of the output file
*/
template<typename RealType, typename TableType>
void write_csv_1d_table(const TableType& table,
    const RealType left, const RealType right,
    const int how_many,
    bool log_scale, const std::string& file_name)
{
    auto coords = std::vector<RealType>(how_many);
    if(log_scale){
            std::generate(coords.begin(), coords.end(), [=,i = 0]() mutable{
            return std::exp(
                std::log(left) + (i++)*(std::log(right)-std::log(left))/(how_many-1));
        });
    }
    else
    {
        std::generate(coords.begin(), coords.end(), [=,i = 0]() mutable{
            return left + (i++)*(right-left)/(how_many-1);
        });
    }

    auto res = std::vector<RealType>(how_many);
#ifdef PXRMP_TABLE_GEN_HAS_OPENMP
    #pragma omp parallel
#endif
    for(int i = 0 ; i < how_many; ++i){
        res[i] = table.interp(coords[i]);
    }

    std::ofstream of{file_name};
    for (int i = 0; i < how_many; ++i){
        of << coords[i] << ", " <<  res[i] << "\n";
    }
    of.close();
}

/**
* Writes a 2D lookup table in csv format
*
* @tparam RealType the floating point type to be used
* @tparam TableType the table type to be used
* @param[in] x1 the left limit of the first coordinate to sample the table
* @param[in] x2 the right limit of the first coordinate to sample the table
* @param[in] y1 the left limit of the second coordinate to sample the table
* @param[in] y2 the right limit of the second coordinate to sample the table
* @param[in] how_many_x controls how many points are used for the first coordinate to sample the lookup table
* @param[in] how_many_y controls how many points are used for the second coordinate to sample the lookup table
* @param[in] log_scale_x controls if the first axis of the table is in log or linear scale
* @param[in] log_scale_y controls if the second axis of the table is in log or linear scale
* @param[in] file_name the name of the output file
*/
template<typename RealType, typename TableType>
void write_csv_2d_table(const TableType& table,
    const RealType x1, const RealType x2,
    const RealType y1, const RealType y2,
    const int how_many_x, const int how_many_y,
    bool log_scale_x, bool log_scale_y, const std::string& file_name)
{
    auto coords_x = std::vector<RealType>(how_many_x);
    auto coords_y = std::vector<RealType>(how_many_y);
    if(log_scale_x){
            std::generate(coords_x.begin(), coords_x.end(), [=,i = 0]() mutable{
            return std::exp(
                std::log(x1) + (i++)*(std::log(x2)-std::log(x1))/(how_many_x-1));
        });
    }
    else
    {
        std::generate(coords_x.begin(), coords_x.end(), [=,i = 0]() mutable{
            return x1 + (i++)*(x2-x1)/(how_many_x-1);
        });
    }

    if(log_scale_y){
            std::generate(coords_y.begin(), coords_y.end(), [=,i = 0]() mutable{
            return std::exp(
                std::log(y1) + (i++)*(std::log(y2)-std::log(y1))/(how_many_y-1));
        });
    }
    else
    {
        std::generate(coords_y.begin(), coords_y.end(), [=,i = 0]() mutable{
            return y1 + (i++)*(y2-y1)/(how_many_y-1);
        });
    }

    auto res = std::vector<RealType>(how_many_x * how_many_y);
#ifdef PXRMP_TABLE_GEN_HAS_OPENMP
    #pragma omp parallel
#endif
    for(int i = 0 ; i < how_many_x; ++i){
        for(int j = 0 ; j < how_many_y; ++j){
            res[i*how_many_y + j] = table.interp(coords_x[i], coords_y[j]);
        }
    }

    std::ofstream of{file_name};
    for(int i = 0 ; i < how_many_x; ++i){
        for(int j = 0 ; j < how_many_y; ++j){
            of << coords_x[i] << ", " <<  coords_y[j]  << ", " << res[i*how_many_y+j]/coords_x[i] << "\n";
        }
    }
    of.close();
}

/**
* Generates a dNdt lookup table for the Breit-Wheeler process
*
* @tparam RealType the floating point type to be used
* @tparam Policy the table generation policy (it can be used to force calculations in double precision with a single precision table)
* @param[in] bw_params table parameters
* @param[in] file_name_prefix the prefix of the output file
*/
template<
    typename RealType,
    px_bw::generation_policy Policy = px_bw::generation_policy::regular>
void generate_breit_wheeler_dndt_table(
    px_bw::dndt_lookup_table_params<RealType> bw_params,
    const std::string& file_name_prefix)
{
    auto table = px_bw::dndt_lookup_table<
        RealType, std::vector<RealType>>{bw_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name_prefix + "_dndt.bin"};
    of.write (raw_data.data(),raw_data.size());
    of.close();

    write_csv_1d_table(table, bw_params.chi_phot_min*0.1f, bw_params.chi_phot_max*10.0f,
        bw_params.chi_phot_how_many*10, true, file_name_prefix + "_dndt.csv");
}

/**
* Generates a pair production lookup table for the Breit-Wheeler process
*
* @tparam RealType the floating point type to be used
* @tparam Policy the table generation policy (it can be used to force calculations in double precision with a single precision table)
* @param[in] bw_params table parameters
* @param[in] file_name_prefix the prefix of the output file
*/
template<
    typename RealType,
    px_bw::generation_policy Policy = px_bw::generation_policy::regular>
void generate_breit_wheeler_pair_prod_table(
    px_bw::pair_prod_lookup_table_params<RealType> bw_params,
    const std::string& file_name)
{
    auto table = px_bw::pair_prod_lookup_table<
        RealType, std::vector<RealType>>{
            bw_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name + "_pairprod.bin"};
    of.write (raw_data.data(),raw_data.size());
    of.close();


    write_csv_2d_table(table, bw_params.chi_phot_min*0.1f, bw_params.chi_phot_max*10.f,
        RealType(0.0), RealType(1.0)-std::numeric_limits<RealType>::epsilon(), bw_params.chi_phot_how_many*3,
        bw_params.frac_how_many*3, true, false, file_name + "_pairprod.csv");

}

/**
* Generates a dNdt lookup table for the Quantum Synchrotron process
*
* @tparam RealType the floating point type to be used
* @tparam Policy the table generation policy (it can be used to force calculations in double precision with a single precision table)
* @param[in] qs_params table parameters
* @param[in] file_name_prefix the prefix of the output file
*/
template<
    typename RealType,
    px_qs::generation_policy Policy = px_qs::generation_policy::regular>
void generate_quantum_sync_dndt_table(
     px_qs::dndt_lookup_table_params<RealType> qs_params,
    const std::string& file_name)
{
    auto table = px_qs::dndt_lookup_table<
        RealType, std::vector<RealType>>{
            qs_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name + "_dndt.bin"};
    of.write (raw_data.data(),raw_data.size());
    of.close();

    write_csv_1d_table(table, qs_params.chi_part_min*0.1f, qs_params.chi_part_max*10.0f,
        qs_params.chi_part_how_many*10, true, file_name + "_dndt.csv");
}

/**
* Generates a photon emission lookup table for the Quantum Synchrotron process
*
* @tparam RealType the floating point type to be used
* @tparam Policy the table generation policy (it can be used to force calculations in double precision with a single precision table)
* @param[in] qs_params table parameters
* @param[in] file_name_prefix the prefix of the output file
*/
template<
    typename RealType,
    px_qs::generation_policy Policy = px_qs::generation_policy::regular>
void generate_quantum_sync_photem_table(
    px_qs::photon_emission_lookup_table_params<RealType> qs_params,
    const std::string& file_name)
{
    auto table = px_qs::photon_emission_lookup_table<
        RealType, std::vector<RealType>>{
            qs_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name + "_photem.bin"};
    of.write (raw_data.data(),raw_data.size());
    of.close();

    write_csv_2d_table(table, qs_params.chi_part_min*0.1f, qs_params.chi_part_max*10.f,
        std::numeric_limits<RealType>::epsilon(), RealType(1.0), qs_params.chi_part_how_many*3,
        qs_params.frac_how_many*3, true, false, file_name + "_photem.csv");
}

/**
* Wraps the call to std::stod, checking if the string has been correctly parsed.
* If conversion fails, it prints an error message and exits the program.
*
* @param[in] str the string to be parsed
* @return std::stod(str) if conversion succeeds
*/
double stod_wrapper(const std::string& str)
{
    size_t idx = 0;
    double val = std::stod(str, &idx);
    if (idx != str.length()){
        print_error("Failed to parse '" + str + "' as a floating point number!");
        exit(EXIT_FAILURE);
    }
    return val;
}

/**
* Wraps the call to std::stoi, checking if the string has been correctly parsed.
* If conversion fails, it prints an error message and exits the program.
*
* @param[in] str the string to be parsed
* @return std::stoi(str) if conversion succeeds
*/
int stoi_wrapper(const std::string& str)
{
    size_t idx = 0;
    int val = std::stoi(str, &idx);
    if (idx != str.length()){
        print_error("Failed to parse '" + str + "' as an integer number!");
        exit(EXIT_FAILURE);
    }
    return val;
}

/**
* Parses arguments to set parameters for Breit-Wheeler lookup tables generation
*
* @tparam RealType the floating point type to be used
* @param[in] args the command line arguments
* @return a struct containing the parameters to generate Breit-Wheeler lookup tables
*/
template<typename RealType>
BreitWheelerTableParams<RealType>
parse_breit_wheeler_params(std::map<std::string, std::string> args)
{
    BreitWheelerTableParams<RealType> params;

    params.chi_min = px_bw::default_chi_phot_min<RealType>;
    params.chi_max = px_bw::default_chi_phot_max<RealType>;
    params.chi_size = px_bw::default_chi_phot_how_many;
    params.frac_size = px_bw::default_frac_how_many;

    auto s_cmin = args.find(CMD_CHI_MIN);
    if(s_cmin != args.end())
        params.chi_min = static_cast<RealType>(stod_wrapper(s_cmin->second));

    auto s_cmax = args.find(CMD_CHI_MAX);
    if( s_cmax != args.end())
        params.chi_max = static_cast<RealType>(stod_wrapper(s_cmax->second));

    auto s_csize = args.find(CMD_CHI_SIZE);
    if( s_csize != args.end())
        params.chi_size = stoi_wrapper(s_csize->second);

    auto s_fsize = args.find(CMD_FRAC_SIZE);
    if(s_fsize != args.end())
        params.frac_size = stoi_wrapper(s_fsize->second);

    return params;
}

/**
* Wraps the call to Breit-Wheeler lookup table generators and prints some useful pieces
* of information on screen
*
* @tparam RealType the floating point type to be used
* @tparam Policy the table generation policy (it can be used to force calculations in double precision with a single precision table)
* @param[in] params the parameters to be used to generate the tables
* @param[in] file_name_prefix the prefix of the output file
*/
template<
    typename RealType,
    px_bw::generation_policy Policy>
void do_breit_wheeler(BreitWheelerTableParams<RealType> params, const std::string& file_name_prefix)
{
    std::cout << " ***** QED table generator: Breit-Wheeler tables *****" << std::endl;
#ifdef PXRMP_TABLE_GEN_HAS_OPENMP
    std::cout << " Using " << omp_get_max_threads() << " OpenMP threads." << std::endl;
#else
    std::cout << " No OpenMP support." << std::endl;
#endif

    if (std::is_same<RealType, double>::value){
        std::cout << " Tables will be generated in double precision." << std::endl;
    }
    else{
        const auto prec =
            (Policy == px_bw::generation_policy::force_internal_double) ?
            "double" : "single";
        std::cout << " Tables will be calculated in " << prec <<
            " precision and will be saved in single precision." << std::endl;
    }

    std::cout << " Table parameters:\n"
        << "  - chi min   : " << params.chi_min << "\n"
        << "  - chi_max   : " << params.chi_max << "\n"
        << "  - chi_size  : " << params.chi_size << "\n"
        << "  - frac_size : " << params.frac_size << "\n";

    const auto dndt_params =
        px_bw::dndt_lookup_table_params<RealType>{
            params.chi_min,
            params.chi_max,
            params.chi_size};

    const auto pair_prod_params =
        px_bw::pair_prod_lookup_table_params<RealType>{
            params.chi_min,
            params.chi_max,
            params.chi_size,
            params.frac_size};

    generate_breit_wheeler_dndt_table<RealType, Policy>(dndt_params, file_name_prefix);
    generate_breit_wheeler_pair_prod_table<RealType, Policy>(pair_prod_params, file_name_prefix);
}

/**
* Parses arguments to set parameters for Quantum Synchrotron lookup tables generation
*
* @tparam RealType the floating point type to be used
* @param[in] args the command line arguments
* @return a struct containing the parameters to generate Quantum Synchrotron lookup tables
*/
template<typename RealType>
QuantumSyncTableParams<RealType>
parse_quantum_sync_params(std::map<std::string, std::string> args)
{
    QuantumSyncTableParams<RealType> params;

    params.chi_min = px_qs::default_chi_part_min<RealType>;
    params.chi_max = px_qs::default_chi_part_max<RealType>;
    params.chi_size = px_qs::default_chi_part_how_many;
    params.frac_min = px_qs::default_frac_min<RealType>;
    params.frac_size = px_qs::default_frac_how_many;

    auto s_cmin = args.find(CMD_CHI_MIN);
    if(s_cmin != args.end())
        params.chi_min = static_cast<RealType>(stod_wrapper(s_cmin->second));

    auto s_cmax = args.find(CMD_CHI_MAX);
    if(s_cmax != args.end())
        params.chi_max = static_cast<RealType>(stod_wrapper(s_cmax->second));

    auto s_fmin = args.find(CMD_FRAC_MIN);
    if(s_fmin != args.end())
        params.frac_min = static_cast<RealType>(stod_wrapper(s_fmin->second));

    auto s_csize = args.find(CMD_CHI_SIZE);
    if(s_csize != args.end())
        params.chi_size = stoi_wrapper(s_csize->second);

    auto s_fsize = args.find(CMD_FRAC_SIZE);
    if(s_fsize != args.end())
        params.frac_size = stoi_wrapper(s_fsize->second);

    return params;
}

/**
* Wraps the call to Quantum Synchrotron lookup table generators and prints some useful pieces
* of information on screen
*
* @tparam RealType the floating point type to be used
* @tparam Policy the table generation policy (it can be used to force calculations in double precision with a single precision table)
* @param[in] params the parameters to be used to generate the tables
* @param[in] file_name_prefix the prefix of the output file
*/
template<
    typename RealType,
    px_qs::generation_policy Policy>
void do_quantum_sync(QuantumSyncTableParams<RealType> params, const std::string& file_name_prefix)
{
    std::cout << " ***** QED table generator: Quantum Synchrotron tables *****" << std::endl;
#ifdef PXRMP_TABLE_GEN_HAS_OPENMP
    std::cout << " Using " << omp_get_max_threads() << " OpenMP threads." << std::endl;
#else
    std::cout << " No OpenMP support." << std::endl;
#endif

    if (std::is_same<RealType, double>::value){
        std::cout << " Tables will be generated in double precision." << std::endl;
    }
    else{
        const auto prec =
            (Policy == px_qs::generation_policy::force_internal_double) ?
            "double" : "single";
        std::cout << " Tables will be calculated in " << prec <<
            " precision and will be saved in single precision." << std::endl;
    }

    std::cout << " Table parameters:\n"
        << "  - chi min   : " << params.chi_min << "\n"
        << "  - chi_max   : " << params.chi_max << "\n"
        << "  - frac_min  : " << params.frac_min << "\n"
        << "  - chi_size  : " << params.chi_size << "\n"
        << "  - frac_size : " << params.frac_size << "\n";

    const auto dndt_params =
        px_qs::dndt_lookup_table_params<RealType>{
            params.chi_min,
            params.chi_max,
            params.chi_size};

    const auto pair_prod_params =
        px_qs::photon_emission_lookup_table_params<RealType>{
            params.chi_min,
            params.chi_max,
            params.frac_min,
            params.chi_size,
            params.frac_size};

    generate_quantum_sync_dndt_table<RealType, Policy>(dndt_params, file_name_prefix);
    generate_quantum_sync_photem_table<RealType, Policy>(pair_prod_params, file_name_prefix);

    std::cout << " ____________________________" << std::endl;
}

/**
* Prints default values of optional parameters in the help message
*/
void print_default_values()
{  
    std::cout << "   Breit-Wheeler tables: " << std::endl;
    std::cout << "       " << CMD_CHI_MIN << " : "
        << px_bw::default_chi_phot_min<double> << std::endl;
    std::cout << "       " << CMD_CHI_MAX << " : "
        << px_bw::default_chi_phot_max<double> << std::endl;
    std::cout << "       " << CMD_FRAC_MIN << " : "
        << "(unused)" << std::endl;
    std::cout << "       " << CMD_CHI_SIZE << " : "
        << px_bw::default_chi_phot_how_many << std::endl;
    std::cout << "       " << CMD_FRAC_SIZE << " : "
        << px_bw::default_frac_how_many << std::endl;

    std::cout << "   Quantum Synchrotron tables: " << std::endl;
    std::cout << "       " << CMD_CHI_MIN << " : "
        << px_qs::default_chi_part_min<double> << std::endl;
    std::cout << "       " << CMD_CHI_MAX << " : "
        << px_qs::default_chi_part_max<double> << std::endl;
    std::cout << "       " << CMD_FRAC_MIN << " : "
        << px_qs::default_frac_min<double> << std::endl;
    std::cout << "       " << CMD_CHI_SIZE << " : "
        << px_qs::default_chi_part_how_many << std::endl;
    std::cout << "       " << CMD_FRAC_SIZE << " : "
        << px_qs::default_frac_how_many << std::endl;
}

/**
* Prints a help message
*/
void print_help_message()
{
    // Max size of a CMD field in the help message
    const int MAX_CMD_SIZE = 26;

    std::cout << " ***** QED table generator HELP *****" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_HELP_S <<
        " : prints this help message" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_HELP_L <<
        " : prints this help message" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_TABLE <<
        " : sets table type. Must be either " << OPT_TABLE_BREIT_WHEELER <<
        " or " << OPT_TABLE_QUANTUM_SYNC << "." << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_PRECISION <<
        " : sets table precision. Must be one of the following: "
        << OPT_PRECISION_DOUBLE << ", "
        << OPT_PRECISION_SINGLE << ", "
        << OPT_PRECISION_SINGLE_WITH_DOUBLE_COMP << ".\n"
        << std::setw(MAX_CMD_SIZE+3) << "" << OPT_PRECISION_SINGLE_WITH_DOUBLE_COMP
        << " means that calculations are carried out in double precision\n"
        << std::setw(MAX_CMD_SIZE+3) << ""
        << "even if the table is in single precision. " << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_CHI_MIN <<
        " : sets minimum chi parameter (real number, optional*)" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_CHI_MAX <<
        " : sets maximum chi parameter (real number, optional*)" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_FRAC_MIN <<
        " : sets minimum frac parameter (real number, optional*)" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_CHI_SIZE <<
        " : sets number of chi parameters (integer number, optional*)" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_FRAC_SIZE <<
        " : sets number of frac parameters (integer number, optional*)" << std::endl;

    std::cout << std::setw(MAX_CMD_SIZE) << CMD_FILENAME <<
        " : sets output file name (string)" << std::endl;

    std::cout << " * Optional parameters have default values: " << std::endl;
    print_default_values();

    std::cout << " ____________________________" << std::endl;
}

/**
* Checks the number of command line arguments and exits if necessary
*
* @param[in] argc the number of command line arguments
*/
void check_argc(int argc)
{
    // We must pass an even number of arguments (and at least 1 argument),
    // thus argc must be odd and > 2
    if (argc%2 != 1 || argc < 2){
        print_error("The number of arguments is incorrect!");
        print_help_message();
        exit(EXIT_FAILURE);
    }
}

/**
* Handles the special case in which the user asks to print the help message
*
* @param[in] argc the number of command line arguments
* @param[in] argv the command line arguments
*/
void handle_help_argument(int argc, char** argv)
{
    if(argc < 2)
        return;

    const auto arg = std::string{argv[1]};
    if(arg.compare(CMD_HELP_S) == 0 || arg.compare(CMD_HELP_L) == 0){
        print_help_message();
        exit(EXIT_SUCCESS);
    }
}

/**
* Reads command line arguments in (command, argument) pairs
*
* @param[in] argc the number of command line arguments
* @param[in] argv the command line arguments
* @return the command line arguments in a (command --> argument) map
*/
std::map<std::string, std::string> read_arg_pairs (int argc, char** argv)
{
    check_argc(argc);

    auto arg_map = std::map<std::string, std::string>{};

    for (int i = 1; i < argc; i += 2){
        if ( arg_map.find(argv[i]) != arg_map.end() ){
            print_error("Command line argument '" + 
                std::string{argv[i]} + "' is duplicated!" );
            exit(EXIT_FAILURE);
        }
        arg_map[argv[i]] = argv[i+1];
    }

    return arg_map;
}

/**
* Parses the command line arguments to determine table precision
*
* @param[in] args command line arguments organised as a (command --> argument) map
* @return table precision
*/
Precision parse_precision(std::map<std::string, std::string>& args)
{
    Precision precision;
    if(args[CMD_PRECISION] == OPT_PRECISION_DOUBLE){
        precision = Precision::double_precision;
    }
    else if(args[CMD_PRECISION] == OPT_PRECISION_SINGLE){
        precision = Precision::single_precision;
    }
    else if(args[CMD_PRECISION] == OPT_PRECISION_SINGLE_WITH_DOUBLE_COMP){
        precision = Precision::single_prec_out_double_prec_comp;
    }
    else{
        print_error("Precision is wrong or not set!");
        print_help_message();
        exit(EXIT_FAILURE);
    }

    return precision;
}

/**
* Parses the command line arguments to determine table type
*
* @param[in] args command line arguments organised as a (command --> argument) map
* @return table type
*/
TableType parse_table_type(std::map<std::string, std::string>& args){
    TableType table_type;
    if(args[CMD_TABLE] == OPT_TABLE_BREIT_WHEELER){
        table_type = TableType::breit_wheeler_table;
    }
    else if(args[CMD_TABLE] == OPT_TABLE_QUANTUM_SYNC){
        table_type = TableType::quantum_synchrotron_table;
    }
    else{
        print_error("Lookup table type not recognized!");
        print_help_message();
        exit(EXIT_FAILURE);
    }

    return table_type;
}

/**
* Parses the command line arguments to determine output name prefix
*
* @param[in] args command line arguments organised as a (command --> argument) map
* @return output name prefix
*/
std::string parse_file_name_prefix(std::map<std::string, std::string>& args){
    auto fname = args.find(CMD_FILENAME);
    if(fname != args.end()){
        return fname->second;
    }
    else{
        print_error("Filename not provided!");
        print_help_message();
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char** argv)
{
    handle_help_argument(argc, argv);
    auto args = read_arg_pairs(argc, argv);

    const auto precision = parse_precision(args);
    const auto table_type = parse_table_type(args);
    const auto file_name_prefix = parse_file_name_prefix(args);

    if(table_type == TableType::breit_wheeler_table){
        if(precision == Precision::double_precision){
            auto params = parse_breit_wheeler_params<double>(args);
            do_breit_wheeler<
                double, px_bw::generation_policy::regular>(
                    params, file_name_prefix);
        }
        else if(precision == Precision::single_precision){
            auto params = parse_breit_wheeler_params<float>(args);
            do_breit_wheeler<
                float, px_bw::generation_policy::regular>(
                    params, file_name_prefix);
        }
        else if(precision == Precision::single_prec_out_double_prec_comp){
            auto params = parse_breit_wheeler_params<float>(args);
            do_breit_wheeler<
                float, px_bw::generation_policy::force_internal_double>(
                    params, file_name_prefix);
        }
    }
    else if(table_type == TableType::quantum_synchrotron_table){
        if(precision == Precision::double_precision){
            auto params = parse_quantum_sync_params<double>(args);
            do_quantum_sync<
                double, px_qs::generation_policy::regular>(
                    params, file_name_prefix);
        }
        else if(precision == Precision::single_precision){
            auto params = parse_quantum_sync_params<float>(args);
            do_quantum_sync<
                float, px_qs::generation_policy::regular>(
                    params, file_name_prefix);
        }
        else if(precision == Precision::single_prec_out_double_prec_comp){
            auto params = parse_quantum_sync_params<float>(args);
            do_quantum_sync<
                float, px_qs::generation_policy::force_internal_double>(
                    params, file_name_prefix);
        }
    }

    exit(EXIT_SUCCESS);
}
