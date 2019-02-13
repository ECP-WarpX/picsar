#ifndef __PHOTONS__
#define __PHOTONS__

#include <string>
#include <cmath>

#include "commons.h"
#include "species.h"

namespace testbed{
  class photons : public species{
  public:
    photons(std::string name);
    ~photons();
    void push_positions(picsar::ttime dt);
    void push_momenta(picsar::ttime dt);

    //friend photons operator+(const photons &p1, const photons &cp2);

  };
}


#endif
