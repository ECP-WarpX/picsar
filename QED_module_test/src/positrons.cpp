#include "positrons.h"

using namespace testbed;
using namespace picsar;


positrons::positrons(std::string name):leptons(name){
    header = std::string("positrons");
    init_mass_and_charge();
}

positrons::~positrons(){}

void positrons::init_mass_and_charge(){
    mass = 1.0;
    charge = 1.0;
}
