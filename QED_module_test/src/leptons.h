#ifndef __LEPTONS__
#define __LEPTONS__

#include <cmath>

#include "species.h"
#include "commons.h"

namespace testbed{
  class leptons : public species{
  public:
    leptons(std::string name);
    ~leptons();
    void push_positions(time dt);
    void push_momenta(time dt);
  protected:
    virtual void init_mass_and_charge() = 0;
  };
}


#endif
