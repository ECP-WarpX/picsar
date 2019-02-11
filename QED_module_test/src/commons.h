#ifndef __COMMONS__
#define __COMMONS__

#include <vector>
#include <array>
#include <functional>

namespace testbed{
  typedef double time;
  typedef std::array<double, 3> position;
  typedef std::array<double, 3> momentum;
  typedef std::array<std::vector<double>, 3> positions_list;
  typedef std::array<std::vector<double>, 3> momenta_list;

  typedef std::array<double, 6> em_field;
  typedef std::function<em_field(position, time)> em_field_function;
  typedef std::array<std::vector<double>, 6> em_field_list;
}






#endif
