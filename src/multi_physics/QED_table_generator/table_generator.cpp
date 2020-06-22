
#include <string>
#include <vector>
#include <fstream>

/*#include <mpi.h>*/
#include <omp.h>

#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tabulated_functions.hpp"

namespace px_bw = picsar::multi_physics::phys::breit_wheeler;

template<typename RealType>
void generate_breit_wheeler_dndt_table_log(const std::string& file_name)
{
    px_bw::dndt_lookup_table_params<RealType> bw_params{1e-3,1e3,100};

    auto table = px_bw::dndt_lookup_table<
        RealType, std::vector<RealType>,
        px_bw::dndt_table_type::log, px_bw::dndt_table_out_policy::approx>{
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

int main(int argc, char** argv)
{
    //MPI_Init(argc,argv);

    auto file_name = "bw_dndt.tab";
    generate_breit_wheeler_dndt_table_log<float>(file_name);

    //MPI_Finalize();
    return 0;
}
