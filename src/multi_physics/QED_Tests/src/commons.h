#ifndef __COMMONS__
#define __COMMONS__

#include <vector>
#include <array>
#include <functional>

namespace testbed{
  typedef double ttime;
  typedef std::array<double, 3> position;
  typedef std::array<double, 3> momentum;
  typedef std::array<std::vector<double>, 3> positions_list;
  typedef std::array<std::vector<double>, 3> momenta_list;

  typedef std::array<double, 6> em_field;
  typedef std::function<em_field(position, ttime)> em_field_function;
  typedef std::array<std::vector<double>, 6> em_field_list;

  typedef std::function<void(momenta_list&, const em_field_list&, double mass, double charge, ttime dt)> mom_pusher_function;

  typedef std::function<void(positions_list&, momenta_list&,
    const em_field_list&, std::vector<double>&, double mass, double charge, ttime dt)> simple_process;
}

#endif
