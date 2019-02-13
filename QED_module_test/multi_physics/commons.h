#ifndef __PMP_COMMONS__
#define __PMP_COMMONS__

#include <vector>
#include <array>
#include <functional>

namespace picsar{
  typedef double ttime;
  typedef std::array<double, 3> position;
  typedef std::array<double, 3> momentum;
  typedef std::array<std::vector<double>, 3> positions_list;
  typedef std::array<std::vector<double>, 3> momenta_list;

  typedef std::array<double, 6> em_field;
  typedef std::function<em_field(position, ttime)> em_field_function;
  typedef std::array<std::vector<double>, 6> em_field_list;

  typedef std::function<void(momenta_list&, const em_field_list&, double mass, double charge, ttime)> mom_pusher_function;

  const double _km = 1.0e3;
  const double _m = 1.0;
  const double _mm = 1.0e-3;
  const double _um = 1.0e-6;
  const double _nm = 1.0e-9;

}






#endif
