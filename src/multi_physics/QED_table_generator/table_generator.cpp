
#include <string>
#include <vector>
#include <fstream>
#include <chrono>

#include <omp.h>

#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tables.hpp"
#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tabulated_functions.hpp"
#include "../QED/src/utils/progress_bar.hpp"

namespace px_bw = picsar::multi_physics::phys::breit_wheeler;
namespace px_qs = picsar::multi_physics::phys::quantum_sync;
namespace px_ut = picsar::multi_physics::utils;

template<typename RealType>
void generate_breit_wheeler_dndt_table(const std::string& file_name)
{
    px_bw::dndt_lookup_table_params<RealType> bw_params{1e-4,1e4,256};

    auto table = px_bw::dndt_lookup_table<
        RealType, std::vector<RealType>>{bw_params};

    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());

    auto t_start =  std::chrono::system_clock::now();
    int count = 0;
    #pragma omp parallel for
    for (int i = 0; i < all_vals.size(); ++i){
        all_vals[i] = px_bw::compute_T_function(all_coords[i]);
        #pragma omp critical
        {
            count++;
            px_ut::draw_progress(count, all_vals.size(), "BW dndt", 10);
        }
    }
    auto t_end =  std::chrono::system_clock::now();
    px_ut::draw_progress(count, count, "BW dndt", 10, true);
    std::cout << "Done in " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(
            t_end - t_start).count()/1000.0 << " seconds! \n" << std::endl;

    if(!table.set_all_vals(all_vals) )
        throw "Fail!";

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}

template<typename RealType>
void generate_breit_wheeler_pair_prod_table(const std::string& file_name)
{
    const int chi_size = 256;
    const int frac_size = 256;
    px_bw::pair_prod_lookup_table_params<RealType> bw_params{1e-4,1e4,chi_size,frac_size};

    auto table = px_bw::pair_prod_lookup_table_logchi_linfrac<
        RealType, std::vector<RealType>>{
            bw_params};

    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());

    auto fracs = std::vector<RealType>(frac_size);
    for(int j = 0; j < frac_size; ++j)
        fracs[j] = all_coords[j][1];

    auto t_start =  std::chrono::system_clock::now();
    int count = 0;
    #pragma omp parallel for
    for (int i = 0; i < chi_size; ++i){
        const auto temp = px_bw::compute_cumulative_prob(
            all_coords[i*frac_size][0],fracs);

        std::copy(temp.begin(), temp.end(), all_vals.begin()+i*frac_size);

        #pragma omp critical
        {
            count++;
            px_ut::draw_progress(count, chi_size, "BW pair prod", 10);
        }
    }
    auto t_end =  std::chrono::system_clock::now();
    px_ut::draw_progress(count, count, "BW pair prod", 10, true);
    std::cout << "Done in " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(
            t_end - t_start).count()/1000.0 << " seconds! \n" << std::endl;

    if(!table.set_all_vals(all_vals) )
        throw "Fail!";

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}

template<typename RealType>
void generate_quantum_sync_dndt_table(const std::string& file_name)
{
    px_qs::dndt_lookup_table_params<RealType> qs_params{1e-4,1e4,256};

    auto table = px_qs::dndt_lookup_table<
        RealType, std::vector<RealType>,
        px_qs::dndt_table_out_policy::approx>{
            qs_params};

    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());

    auto t_start =  std::chrono::system_clock::now();
    int count = 0;
    #pragma omp parallel for
    for (int i = 0; i < all_vals.size(); ++i){
        all_vals[i] = px_qs::compute_G_function(all_coords[i]);

        #pragma omp critical
        {
            count++;
            px_ut::draw_progress(count, all_vals.size(), "QS dndt", 10);
        }
    }
    px_ut::draw_progress(count, count, "QS dndt", 10, true);
    auto t_end =  std::chrono::system_clock::now();
    std::cout << "Done in " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(
            t_end - t_start).count()/1000.0 << " seconds! \n" << std::endl;


    if(!table.set_all_vals(all_vals) )
        throw "Fail!";

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();
}

template<typename RealType>
void generate_quantum_sync_photem_table(const std::string& file_name)
{
    const int chi_size = 256;
    const int frac_size = 256;
    px_qs::photon_emission_lookup_table_params<RealType> qs_params{
        1e-4,1e4,1e-5,chi_size,frac_size};

    auto table = px_qs::photon_emission_lookup_table_logchi_logfrac<
        RealType, std::vector<RealType>>{
            qs_params};

    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());

    auto fracs = std::vector<RealType>(frac_size);
    for(int j = 0; j < frac_size; ++j)
        fracs[j] = all_coords[j][1];

    auto t_start =  std::chrono::system_clock::now();
    int count = 0;
    #pragma omp parallel for
    for (int i = 0; i < chi_size; ++i){
        const auto temp = px_qs::compute_cumulative_prob(
            all_coords[i*frac_size][0],fracs);

        std::copy(temp.begin(), temp.end(), all_vals.begin()+i*frac_size);

        #pragma omp critical
        {
            count++;
            px_ut::draw_progress(count, chi_size, "QS phot em", 10);
        }
    }
    auto t_end =  std::chrono::system_clock::now();
    px_ut::draw_progress(count, chi_size, "QS phot em", 10, true);
    std::cout << "Done in " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(
            t_end - t_start).count()/1000.0 << " seconds! \n" << std::endl;

    if(!table.set_all_vals(all_vals) )
        throw "Fail!";

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}


int main(int argc, char** argv)
{

    std::cout << "** Table generation in double precision ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<double>("bw_dndt.tab");
    generate_breit_wheeler_pair_prod_table<double>("bw_pairprod.tab");
    generate_quantum_sync_dndt_table<double>("qs_dndt.tab");
    generate_quantum_sync_photem_table<double>("qs_photem.tab");

    std::cout << "** Table generation in single precision ** \n" << std::endl;
    generate_breit_wheeler_dndt_table<float>("bw_dndt.ftab");
    generate_breit_wheeler_pair_prod_table<float>("bw_pairprod.ftab");
    generate_quantum_sync_dndt_table<float>("qs_dndt.ftab");
    generate_quantum_sync_photem_table<float>("qs_photem.ftab");

    return 0;
}
