
#include <string>
#include <vector>
#include <fstream>

/*#include <mpi.h>*/
#include <omp.h>

#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tables.hpp"
#include "../QED/src/physics/quantum_sync/quantum_sync_engine_tabulated_functions.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tabulated_functions.hpp"

namespace px_bw = picsar::multi_physics::phys::breit_wheeler;
namespace px_qs = picsar::multi_physics::phys::quantum_sync;

template<typename RealType>
void generate_breit_wheeler_dndt_table_log(const std::string& file_name)
{
    px_bw::dndt_lookup_table_params<RealType> bw_params{1e-4,1e4,256};

    auto table = px_bw::dndt_lookup_table<
        RealType, std::vector<RealType>,
        px_bw::dndt_table_out_policy::approx>{
            bw_params};

    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());

    #pragma omp parallel
    for (int i = 0; i < all_vals.size(); ++i){
        all_vals[i] = px_bw::compute_T_function(all_coords[i]);
    }

    if(!table.set_all_vals(all_vals) )
        throw "Fail!";

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}

template<typename RealType>
void generate_quantum_sync_dndt_table_log(const std::string& file_name)
{
    px_qs::dndt_lookup_table_params<RealType> qs_params{1e-4,1e4,256};

    auto table = px_qs::dndt_lookup_table<
        RealType, std::vector<RealType>,
        px_qs::dndt_table_out_policy::approx>{
            qs_params};

    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());

    #pragma omp parallel
    for (int i = 0; i < all_vals.size(); ++i){
        all_vals[i] = px_qs::compute_G_function(all_coords[i]);
    }

    if(!table.set_all_vals(all_vals) )
        throw "Fail!";

    const auto raw_data = table.serialize();

    std::ofstream of{file_name};
    of.write (raw_data.data(),raw_data.size());
    of.close();

}

int main(int argc, char** argv)
{
    //MPI_Init(argc,argv);

    generate_breit_wheeler_dndt_table_log<double>("bw_dndt.tab");

    generate_quantum_sync_dndt_table_log<double>("qs_dndt.tab");

    //MPI_Finalize();
    return 0;
}
