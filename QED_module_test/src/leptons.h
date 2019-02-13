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
    void push_positions(picsar::ttime dt);
    void push_momenta(picsar::ttime dt);
    void replace_pusher_momenta(picsar::mom_pusher_function pusher);
  protected:
    virtual void init_mass_and_charge() = 0;

    bool is_boris_replaced = false;
    void do_boris_pusher(picsar::ttime dt);
    void do_alternative_push_momenta(picsar::ttime dt);
    picsar::mom_pusher_function pusher;

  };
}


#endif
