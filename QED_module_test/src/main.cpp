#include <iostream>
#include <vector>
#include <memory>
#include <functional>
#include <tuple>
#include <cmath>
#include <cstdint>
#include <fstream>

#include "commons.h"
#include "species.h"
#include "photons.h"
#include "electrons.h"
#include "positrons.h"

#include "landau_lifshitz.h"
#include "nonlin_breit_wheeler_engine.h"
#include "special_functions.h"
#include "quadrature.h"

using namespace std;
using namespace testbed;

const double dt = 0.01;
const int num_steps = 10000;

inline bool is_out(int t_step){return (t_step % 1000 == 0);}

const int64_t seed = 3397169560718639567;

int main(int argc, char** argv){

    cout << "********************Test special functions***************************" << endl;
    cout << "k_1_3 (0.5) = " << picsar::multi_physics::k_1_3 (0.5) << " (exp. 0.989031)" << endl;
    cout << "k_1_3 (1.0) = " << picsar::multi_physics::k_1_3 (1.0) << " (exp. 0.438431)" << endl;
    cout << "k_1_3 (2.0) = " << picsar::multi_physics::k_1_3 (2.0) << " (exp. 0.116545)" << endl;
    cout << "k_2_3 (0.5) = " << picsar::multi_physics::k_2_3 (0.5) << "  (exp. 1.20593 )" << endl;
    cout << "k_2_3 (1.0) = " << picsar::multi_physics::k_2_3 (1.0) << " (exp. 0.494475)" << endl;
    cout << "k_2_3 (2.0) = " << picsar::multi_physics::k_2_3 (2.0) << " (exp. 0.124839)" << endl;
    cout << "*********************************************************************" << endl;
    cout << endl;

    cout << "********************Test numerical integration***********************" << endl;
    cout << "quad(sin(x)**2, 0, 2pi) = " <<
        picsar::multi_physics::quad_a_b([](double x){return sin(x)*sin(x);}, 0, 2*M_PI) << " (exp. 3.14159)" << endl;

    cout << "quad(exp(-x**2), 0, inf.) = " <<
        picsar::multi_physics::quad_a_inf([](double x){return exp(-x*x);},0) << " (exp. 0.886227)" << endl;

    cout << "quad(1/(1+x**2), 0, inf.) = " <<
        picsar::multi_physics::quad_a_inf([](double x){return 1.0/(1+x*x);},0) << " (exp. 1.57079)" << endl;

    cout << "*********************************************************************" << endl;
    cout << endl;

    std::ofstream dump_file{"dump.dat"}; //dump file to store data for later analysis

    cout << "********************QED module testbed***************************" << endl;

    //Fix lambda
    double lambda = 1 * picsar::multi_physics::_um;

    //Init nonlin_breit_wheeler_engine
    picsar::multi_physics::nonlin_breit_wheeler_engine breit_wheeler_engine{seed, lambda};

    vector<shared_ptr<species>> specs;
    //Init a photon
    auto ptr_phot1 = make_shared<photons>("phot1");
    ptr_phot1->add_particle({0,0,0},{0.0,0.0,0.0});
    specs.emplace_back(ptr_phot1);

    //Create LL pusher using multi_physics library
    auto pusher =
    [lambda](momenta_list& mom, const em_field_list& fields,  double mass, double charge, double dt)->void{
        picsar::multi_physics::boris_plus_landau_lifshitz_push(
        mom[0], mom[1], mom[2],
        fields[0], fields[1], fields[2],
        fields[3], fields[4], fields[5],
        mass, charge, dt, lambda);
    };

    //Init an electron
    auto ptr_ele1 = make_shared<electrons>("ele1");
    ptr_ele1->add_particle({0,0,0},{0,100,0.0});
    //Replace pusher
    ptr_ele1->replace_pusher_momenta(pusher);
    specs.emplace_back(ptr_ele1);

    //Init a positron
    auto ptr_pos1 = make_shared<positrons>("pos1");
    ptr_pos1->add_particle({0,0,0},{0.0,0.0,0.0});
    //ptr_pos1->replace_pusher_momenta(pusher);
    specs.emplace_back(ptr_pos1);

    // Main loop
    for (int i = 0; i < num_steps; i++){

    //Dump_file output of a particle position for debug purposes
        auto tpos = ptr_ele1->get_copy_of_positions();
        dump_file << tpos[0][0] << " " << tpos[1][0] << " " << tpos[2][0] << " ";
        auto tmom = ptr_ele1->get_copy_of_momenta();
        dump_file << tmom[0][0] << " " << tmom[1][0] << " " << tmom[2][0] << endl;

        for (auto& sp : specs)
            sp->push_momenta(dt);

        for (auto& sp : specs)
            sp->calc_fields([](position pos, double ttime){return em_field{0,0,0,100.0,0.0,0.0};}, i*dt);

        for (auto& sp : specs)
            sp->push_positions(dt);

        if(is_out(i)){
            for (auto& sp : specs)
                sp->print_on_disk("out", i);
        }
    }

    cout << "*****************************************************************" << endl;

    dump_file.close();

    return 0;
}

//Nice field!
//sp->calc_fields([](position pos, double ttime){return em_field{0,cos(pos[0]),sin(pos[0]),0,sin(pos[0]),-cos(pos[0])};}, i*dt);
