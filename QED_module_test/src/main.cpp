#include <iostream>
#include <vector>
#include <memory>

#include "commons.h"
#include "species.h"
#include "photons.h"
#include "electrons.h"
#include "positrons.h"

using namespace std;
using namespace testbed;

const double dt = 0.001;
const int num_steps = 10000;

bool is_out(int t_step){return (t_step % 1000 == 0);}

int main(int argc, char** argv){
    cout << "********************QED module testbed***************************" << endl;


    vector<shared_ptr<species>> specs;
    //Init a photon
    auto ptr_phot1 = shared_ptr<photons>(new photons("phot1"));
    ptr_phot1->add_particle({0,0,0},{1.0,0.0,0.0});
    specs.emplace_back(ptr_phot1);

    //Init an electron
    auto ptr_ele1 = shared_ptr<electrons>(new electrons("ele1"));
    ptr_ele1->add_particle({0,0,0},{1.0,0.0,0.0});
    specs.emplace_back(ptr_ele1);

    //Init a positron
    auto ptr_pos1 = shared_ptr<positrons>(new positrons("pos1"));
    ptr_pos1->add_particle({0,0,0},{1.0,0.0,0.0});
    specs.emplace_back(ptr_pos1);

    // Main loop
    for (int i = 0; i < num_steps; i++){
        for (auto& sp : specs)
            sp->push_momenta(dt);

        for (auto& sp : specs)
            sp->calc_fields([](position pos, double time){return em_field{0,0,0,0,0,1};}, i*dt);

        for (auto& sp : specs)
            sp->push_positions(dt);

        if(is_out(i)){
            for (auto& sp : specs)
                sp->print_on_disk("out", i);
        }
    }

    cout << "*****************************************************************" << endl;

    return 0;
}
