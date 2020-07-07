
#include <string>
#include <vector>
#include <fstream>
#include <chrono>

#include <omp.h>

#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tables.hpp"
#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"
#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"

namespace px_bw = picsar::multi_physics::phys::breit_wheeler;
namespace px_qs = picsar::multi_physics::phys::quantum_sync;
namespace px_ut = picsar::multi_physics::utils;

template<
    typename RealType,
    px_bw::generation_policy Policy = px_bw::generation_policy::regular>
void generate_breit_wheeler_dndt_table(
    const std::string& file_name)
{
    px_bw::dndt_lookup_table_params<RealType> bw_params{1e-4,1e4,256};

    auto table = px_bw::dndt_lookup_table<
        RealType, std::vector<RealType>>{bw_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}

template<
    typename RealType,
    px_bw::generation_policy Policy = px_bw::generation_policy::regular>
void generate_breit_wheeler_pair_prod_table(
    const std::string& file_name)
{
    const int chi_size = 128;
    const int frac_size = 128;
    px_bw::pair_prod_lookup_table_params<RealType> bw_params{1e-4,1e4,chi_size,frac_size};

    auto table = px_bw::pair_prod_lookup_table_logchi_linfrac<
        RealType, std::vector<RealType>>{
            bw_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}

template<
    typename RealType,
    px_qs::generation_policy Policy = px_qs::generation_policy::regular>
void generate_quantum_sync_dndt_table(
    const std::string& file_name)
{
    px_qs::dndt_lookup_table_params<RealType> qs_params{1e-4,1e4,256};

    auto table = px_qs::dndt_lookup_table<
        RealType, std::vector<RealType>>{
            qs_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();
}

template<
    typename RealType,
    px_qs::generation_policy Policy = px_qs::generation_policy::regular>
void generate_quantum_sync_photem_table(
    const std::string& file_name)
{
    const int chi_size = 128;
    const int frac_size = 128;
    px_qs::photon_emission_lookup_table_params<RealType> qs_params{
        1e-4,1e4,1e-5,chi_size,frac_size};

    auto table = px_qs::photon_emission_lookup_table_logchi_logfrac<
        RealType, std::vector<RealType>>{
            qs_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}


int main(int argc, char** argv)
{

    std::cout << "** Double precision tables ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<double>(
        "bw_dndt.tab");
    generate_breit_wheeler_pair_prod_table<double>(
        "bw_pairprod.tab");
    generate_quantum_sync_dndt_table<double>(
        "qs_dndt.tab");
    generate_quantum_sync_photem_table<double>(
        "qs_photem.tab");

    std::cout << "** Single precision tables calculated in double precision ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<float,
        px_bw::generation_policy::force_internal_double>("bw_dndt.ftab");
    generate_breit_wheeler_pair_prod_table<float,
        px_bw::generation_policy::force_internal_double>("bw_pairprod.ftab");
    generate_quantum_sync_dndt_table<float,
        px_qs::generation_policy::force_internal_double>("qs_dndt.ftab");
    generate_quantum_sync_photem_table<float,
        px_qs::generation_policy::force_internal_double>("qs_photem.ftab");

    std::cout << "** Single precision tables ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<float>(
        "bw_dndt.ftab");
    generate_breit_wheeler_pair_prod_table<float>(
        "bw_pairprod.ftab");
    generate_quantum_sync_dndt_table<float>(
        "qs_dndt.ftab");
    generate_quantum_sync_photem_table<float>(
        "qs_photem.ftab");

    return 0;
}
