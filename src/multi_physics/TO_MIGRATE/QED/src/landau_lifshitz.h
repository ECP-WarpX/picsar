#ifndef __PMP_LANDAU_LIFSHITZ__
#define __PMP_LANDAU_LIFSHITZ__

#include <vector>
#include <cmath>

#include "commons.h"

namespace picsar{
    namespace multi_physics{

        //**********************************************************************
        //This function implements Boris Pusher with Landau Lifshitz algorithm,
        //according to the algorithm proposed by M.Tamburini.
        //It accepts three vectors for the different components of the particle momenta
        //as a refence, 6 vectors for the field components as const references,
        //a mass, a charge, a dt and the reference length lambda.
        //It returns immediately if one of these last 4 parameters is 0.
        //The function gets the number of particles from px vector, thus make sure
        //that all the vectors have the same size (or are at least as big as px).
        //The function assumes that omega = 1 in code units.
        //**********************************************************************
        void boris_plus_landau_lifshitz_push(
        std::vector<double>& px, std::vector<double>& py, std::vector<double>& pz,
        const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez,
        const std::vector<double>& bx, const std::vector<double>& by, const std::vector<double>& bz,
        double mass, double charge, double dt, double lambda);
    }
}

#endif
