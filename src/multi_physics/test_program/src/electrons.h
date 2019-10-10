#ifndef __ELECTRONS__
#define __ELECTRONS__

#include "species.h"
#include "leptons.h"

namespace testbed{
  class electrons : public leptons{
  public:
    electrons(std::string name);
    ~electrons();
  private:
    void init_mass_and_charge();
  };
}



#endif
