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
    void push_positions(ttime dt);
    void push_momenta(ttime dt);
    void replace_pusher_momenta(mom_pusher_function pusher);
  protected:
    virtual void init_mass_and_charge() = 0;

    bool is_boris_replaced = false;
    void do_boris_pusher(ttime dt);
    void do_alternative_push_momenta(ttime dt);
    mom_pusher_function pusher;

  };
}


#endif
