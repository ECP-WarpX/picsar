#ifndef __POSITRONS__
#define __POSITRONS__

#include "species.h"
#include "leptons.h"

namespace testbed{
  class positrons : public leptons{
  public:
    positrons(std::string name);
    ~positrons();
  private:
    void init_mass_and_charge();
  };
}



#endif
