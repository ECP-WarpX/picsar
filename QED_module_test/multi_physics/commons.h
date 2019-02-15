#ifndef __PMP_COMMONS__
#define __PMP_COMMONS__

//#include <vector>
//#include <array>
//#include <functional>

namespace picsar{
    namespace multi_physics{

        //typedef std::function<void(momenta_list&, const em_field_list&, double mass, double charge, ttime)> mom_pusher_function;

        const double _km = 1.0e3;
        const double _m = 1.0;
        const double _cm = 1.0e-2;
        const double _mm = 1.0e-3;
        const double _um = 1.0e-6;
        const double _nm = 1.0e-9;

        const double classical_electron_radius = 2.8179403267e-13 * _cm;
    }
}

#endif
