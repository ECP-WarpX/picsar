//STL
#include <iostream>
#include <vector>
#include <memory>
#include <functional>
#include <tuple>
#include <cmath>
#include <cstdint>
#include <fstream>

//Testbed header files
#include "commons.h"
#include "species.h"
#include "photons.h"
#include "electrons.h"
#include "positrons.h"

//QED library files
#define PXRMP_USE_NORMALIZED_UNITS
#include "../QED/src/breit_wheeler_engine.hpp"


#include <fenv.h>


using namespace std;
using namespace testbed;
namespace pxrmp =  picsar::multi_physics;

//Test cases
void test_BW();

int main(){
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    test_BW();
}

//Helper function
bool does_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}


void test_BW(){
    const double dt = 0.05;
    const int num_steps = 200;

    auto is_out = [](int t_step){return (t_step % 1 == 0);};

    const int64_t seed_BW = 2997169562813639567;

    //Fix lambda to 800 nm
    double lambda = 0.8 * pxrmp::si_micrometer;

    pxrmp::stl_rng_wrapper<double> wrap{seed_BW};

    //This is a struct which contains parameters relevant for BW engine
    pxrmp::breit_wheeler_engine_ctrl<double> bw_ctrl;

    bw_ctrl.chi_phot_tdndt_how_many = 200;

    bw_ctrl.chi_phot_tpair_min = 0.01;
    bw_ctrl.chi_phot_tpair_max = 50;
    bw_ctrl.chi_phot_tpair_how_many = 200;

    bw_ctrl.chi_frac_tpair_how_many = 100;

    //Creates the BW engine
    auto bw_engine =  pxrmp::breit_wheeler_engine<double, pxrmp::stl_rng_wrapper<double>>{
        std::move(wrap),
        lambda, //Optional argument (default value is 1)
        bw_ctrl //Optional argument (with reasonable default values)
        };

    //Generates tables if they do not exist
    if(!does_file_exist("tdndt.bin")){
        bw_engine.compute_dN_dt_lookup_table(&std::cout);
        bw_engine.write_dN_dt_table("tdndt.bin");
    }
    else{
            bw_engine.read_dN_dt_table("tdndt.bin");
    }

    if(!does_file_exist("tpair.bin")){
        bw_engine.compute_cumulative_pair_table(&std::cout);
        bw_engine.write_cumulative_pair_table("tpair.bin");
    }
    else{
            bw_engine.read_cumulative_pair_table("tpair.bin");
    }

    //Test
    auto innards = bw_engine.export_innards();
    auto bw_clone = pxrmp::breit_wheeler_engine<double, pxrmp::stl_rng_wrapper<double>>
                    (innards);

    //Init some photons
    auto ptr_phot1 = make_shared<photons>("phot1");
    const size_t how_many_phot = 20000;
    const int64_t seed_photons = 123043949938;

    pxrmp::stl_rng_wrapper<double> wrap_phot{seed_photons};

    for(size_t i = 0; i < how_many_phot; i++){
        // double mom = wrap_phot.unf<double>(8000,12000);
        // double theta = wrap_phot.unf<double>(0,2.0*M_PI);
        // double phi = acos(wrap_phot.unf<double>(-1.0,1.0));
        //
        // ptr_phot1->add_particle({0,0,0},{mom*sin(phi)*cos(theta), mom*sin(phi)*sin(theta) , mom*cos(phi)});
        double momx = wrap_phot.unf(-8000,8000);
        double momy = wrap_phot.unf(-8000,8000);
        double momz = wrap_phot.unf(-8000,8000);

        ptr_phot1->add_particle({0,0,0},{momx, momy, momz});

    }

    for (auto& opt: ptr_phot1->get_ref_of_optical_depth())
        opt = bw_engine.get_optical_depth();

    //Add BW decrease of optical depth to photons
    auto BW_opticaldepth =
    [&bw_clone](positions_list& pos, momenta_list& mom,
      const em_field_list& fields, std::vector<double>& opt_depth, std::vector<bool>& flag, double, double, ttime dt)->void{
          for(size_t i = 0; i < pos[0].size(); i++){
              bool has_event_happend;
              double dt_prod;
              std::tie(has_event_happend, dt_prod) =
                  bw_clone.evolve_opt_depth_and_determine_event
                  (mom[0][i], mom[1][i], mom[2][i],
                  fields[0][i], fields[1][i], fields[2][i],
                  fields[3][i], fields[4][i], fields[5][i], dt, opt_depth[i]);
              flag[i] = has_event_happend;
              if(has_event_happend){
                  double inv_ptotal =
                    1.0/sqrt(mom[0][i]*mom[0][i]
                    + mom[1][i]*mom[1][i] + mom[2][i]*mom[2][i]);
                  double coeff = dt_prod*inv_ptotal;
                  pos[0][i] += coeff*mom[0][i];
                  pos[1][i] += coeff*mom[1][i];
                  pos[2][i] += coeff*mom[2][i];
              }
          }

    };
    ptr_phot1->add_simple_process(BW_opticaldepth, 0);
    //
    //
        //Add simple process to print on disk how many photons are left
        auto BW_howmanyphotons =
        [&bw_clone, lambda](positions_list& pos, momenta_list& ,
          const em_field_list& , std::vector<double>&, std::vector<bool>& , double, double, ttime dt)->void{
              std::cout << "There are " << pos[0].size() << " photons left!" << std::endl;
        };
        ptr_phot1->add_simple_process(BW_howmanyphotons, 777);

    //
    //Create EMPTY electron and positron species
    auto ptr_ele1 = make_shared<electrons>("ele1");
    auto ptr_pos1 = make_shared<positrons>("pos1");
    //
    //add_process_with_destruction
    auto BW_pair_prod=
    [&bw_clone, lambda, ptr_ele1, ptr_pos1](const position& pos, const momentum& mom,
    const em_field& field, const ooptical_depth& opt, const bool& flag, double mass, double charge, ttime dt)->bool{
        if(!flag)
            return false;

        size_t sampling = 1;
        auto all = bw_clone.generate_breit_wheeler_pairs
            (mom[0], mom[1], mom[2],
            field[0], field[1], field[2],
            field[3], field[4], field[5], 1.0, sampling);

       	auto elec = all[0][0];
       	auto post = all[1][0];

        momentum mom_ele{elec.first[0], elec.first[1], elec.first[2]};
        momentum mom_pos{post.first[0], post.first[1], post.first[2]};
        ptr_ele1->add_particle(pos, mom_ele);
        ptr_pos1->add_particle(pos, mom_pos);

           return true;
       };

    ptr_phot1->add_process_with_destruction(BW_pair_prod, 42);
    //

    //Prepare a species vector
    vector<shared_ptr<species>> specs{ptr_phot1, ptr_ele1, ptr_pos1};

    // Main loop
    for (int i = 0; i < num_steps; i++){
        for (auto& sp : specs)
            sp->calc_fields([](position, double){return em_field{0.0, 0.0, 0.0, 2000.0, 0.0, 0.0};}, i*dt);

        for (auto& sp : specs)
            sp->push_momenta(dt);

        for (auto& sp : specs)
            sp->push_positions(dt);

        for (auto& sp : specs)
            sp->do_all_simple_processes(dt);

        ptr_phot1->do_process_with_destruction(42, dt);

        if(is_out(i)){
            for (auto& sp : specs)
                sp->print_on_disk("out", i);
        }
    }
    cout << "BW END" << endl;
}



//Nice field!
//sp->calc_fields([](position pos, double ttime){return em_field{0,cos(pos[0]),sin(pos[0]),0,sin(pos[0]),-cos(pos[0])};}, i*dt);
