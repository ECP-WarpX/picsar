#include <iostream>
#include <vector>
#include <memory>
#include <functional>

#include "commons.h"
#include "species.h"
#include "photons.h"
#include "electrons.h"
#include "positrons.h"

#include "landau_lifshitz.h"

using namespace std;
using namespace testbed;
using namespace picsar;

const double dt = 0.001;
const int num_steps = 10000;

bool is_out(int t_step){return (t_step % 1000 == 0);}

int main(int argc, char** argv){
    cout << "********************QED module testbed***************************" << endl;

    //Fix lambda
    double lambda = 800 * _nm;

    vector<shared_ptr<species>> specs;
    //Init a photon
    auto ptr_phot1 = make_shared<photons>("phot1");
    ptr_phot1->add_particle({0,0,0},{1.0,0.0,0.0});
    specs.emplace_back(ptr_phot1);

    //Create LL pusher
    auto pusher =
    [lambda](momenta_list& mom, const em_field_list& fields,  double mass, double charge, double dt)->void{
        boris_plus_landau_lifshitz_push(mom, fields, mass, charge, dt, lambda);
        return;
    };

    //Init an electron
    auto ptr_ele1 = make_shared<electrons>("ele1");
    ptr_ele1->add_particle({0,0,0},{1.0,0.0,0.0});
    //Replace pusher
    ptr_ele1->replace_pusher_momenta(pusher);
    specs.emplace_back(ptr_ele1);

    //Init a positron
    auto ptr_pos1 = make_shared<positrons>("pos1");
    ptr_pos1->add_particle({0,0,0},{1.0,0.0,0.0});
    ptr_pos1->replace_pusher_momenta(pusher);
    specs.emplace_back(ptr_pos1);

    // Main loop
    for (int i = 0; i < num_steps; i++){
        for (auto& sp : specs)
            sp->push_momenta(dt);

        for (auto& sp : specs)
            sp->calc_fields([](position pos, double ttime){return em_field{0,0,0,0,0,1};}, i*dt);

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
