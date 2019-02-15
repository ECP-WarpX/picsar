#ifndef __PMP_LANDAU_LIFSHITZ__
#define __PMP_LANDAU_LIFSHITZ__

#include <cmath>

#include "commons.h"

namespace picsar{
    void boris_plus_landau_lifshitz_push(momenta_list& mom, const em_field_list& fields, const double mass, const double charge, const double dt, const double lambda);
}

#endif
