/**
* This program tests the generation of the lookup tables of the QED library and
* shows how to convert them into a binary file.
* For each table it produces also a csv file which can be inspected with pyplot
* or gnuplot.
*/

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <omp.h>

#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tables.hpp"
#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"
#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"

namespace px_bw = picsar::multi_physics::phys::breit_wheeler;
namespace px_qs = picsar::multi_physics::phys::quantum_sync;
namespace px_ut = picsar::multi_physics::utils;

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
    #pragma omp parallel
    for(int i = 0 ; i < how_many; ++i){
        res[i] = table.interp(coords[i]);
    }

    std::ofstream of{file_name};
    for (int i = 0; i < how_many; ++i){
        of << coords[i] << ", " <<  res[i] << "\n";
    }
    of.close();
}

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
    #pragma omp parallel
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

template<
    typename RealType,
    px_bw::generation_policy Policy = px_bw::generation_policy::regular>
void generate_breit_wheeler_dndt_table(
    px_bw::dndt_lookup_table_params<RealType> bw_params,
    const std::string& file_name)
{
    auto table = px_bw::dndt_lookup_table<
        RealType, std::vector<RealType>>{bw_params};

    table.template generate<Policy>();

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

    write_csv_1d_table(table, bw_params.chi_phot_min*0.1f, bw_params.chi_phot_max*10.0f,
        bw_params.chi_phot_how_many*10, true, file_name + ".csv");
}

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

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();


    write_csv_2d_table(table, bw_params.chi_phot_min*0.1f, bw_params.chi_phot_max*10.f,
        RealType(0.0), RealType(1.0)-std::numeric_limits<RealType>::epsilon(), bw_params.chi_phot_how_many*3,
        bw_params.frac_how_many*3, true, false, file_name + ".csv");

}

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

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

    write_csv_1d_table(table, qs_params.chi_part_min*0.1f, qs_params.chi_part_max*10.0f,
        qs_params.chi_part_how_many*10, true, file_name + ".csv");
}

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

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

    write_csv_2d_table(table, qs_params.chi_part_min*0.1f, qs_params.chi_part_max*10.f,
        std::numeric_limits<RealType>::epsilon(), RealType(1.0), qs_params.chi_part_how_many*3,
        qs_params.frac_how_many*3, true, false, file_name + ".csv");
}


int main(int argc, char** argv)
{

    std::cout << "** Double precision tables ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<double>(
        px_bw::default_dndt_lookup_table_params<double>,
        "bw_dndt_d");
    generate_breit_wheeler_pair_prod_table<double>(
        px_bw::default_pair_prod_lookup_table_params<double>,
        "bw_pairprod_d");
    generate_quantum_sync_dndt_table<double>(
        px_qs::default_dndt_lookup_table_params<double>,
        "qs_dndt_d");
    generate_quantum_sync_photem_table<double>(
        px_qs::default_photon_emission_lookup_table_params<double>,
        "qs_photem_d");
    std::cout << "____________________________ \n" << std::endl;

/*
    std::cout << "** Single precision tables calculated in double precision ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<float,
        px_bw::generation_policy::force_internal_double>(
        px_bw::default_dndt_lookup_table_params<float>,
        "bw_dndt_fd");
    generate_breit_wheeler_pair_prod_table<float,
        px_bw::generation_policy::force_internal_double>(
        px_bw::default_pair_prod_lookup_table_params<float>,
        "bw_pairprod_fd");
    generate_quantum_sync_dndt_table<float,
        px_qs::generation_policy::force_internal_double>(
        px_qs::default_dndt_lookup_table_params<float>,
        "qs_dndt_fd");
    generate_quantum_sync_photem_table<float,
        px_qs::generation_policy::force_internal_double>(
        px_qs::default_photon_emission_lookup_table_params<float>,
        "qs_photem_fd");

    std::cout << "____________________________ \n" << std::endl;
    */

    std::cout << "** Single precision tables ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<float>(
        px_bw::default_dndt_lookup_table_params<float>,
        "bw_dndt_f");
    generate_breit_wheeler_pair_prod_table<float>(
        px_bw::default_pair_prod_lookup_table_params<float>,
        "bw_pairprod_f");
    generate_quantum_sync_dndt_table<float>(
        px_qs::default_dndt_lookup_table_params<float>,
        "qs_dndt_f");
    generate_quantum_sync_photem_table<float>(
        px_qs::default_photon_emission_lookup_table_params<float>,
        "qs_photem_f");

    std::cout << "____________________________ \n" << std::endl;

    return 0;
}
