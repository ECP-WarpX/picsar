#include "landau_lifshitz.h"

using namespace picsar;

void picsar::boris_plus_landau_lifshitz_push(momenta_list& mom, const em_field_list& fields, double mass, double charge, double dt, double lambda){
    double beta = 0.5 * charge * dt / mass;

    int num_particles = mom[0].size();
    for (int i = 0; i < num_particles; i++){

        double ptx = mom[0][i] + beta *  fields[0][i];
        double pty = mom[1][i] + beta *  fields[1][i];
        double ptz = mom[2][i] + beta *  fields[2][i];

        double inv_gamman_ttimes_beta = beta/sqrt(mass*mass + mom[0][i]*mom[0][i] + mom[1][i]*mom[1][i] + mom[2][i]*mom[2][i]);

        double bx = inv_gamman_ttimes_beta * fields[3][i];
        double by = inv_gamman_ttimes_beta * fields[4][i];
        double bz = inv_gamman_ttimes_beta * fields[5][i];

        double b2 = bx*bx + by*by + bz*bz;
        double coeff = 1.0 / (1.0 + b2);

        double ptvbx = pty*bz - ptz*by;
        double ptvby = ptz*bx - ptx*bz;
        double ptvbz = ptx*by - pty*bx;

        double ptdb = ptx * bx + pty * by + ptz * bz;

        double bptdbx = ptdb*bx;
        double bptdby = ptdb*by;
        double bptdbz = ptdb*bz;

        double pnx = coeff*(ptx + ptvbx + bptdbx);
        double pny = coeff*(pty + ptvby + bptdby);
        double pnz = coeff*(ptz + ptvbz + bptdbz);

        // ADD BORIS

         mom[0][i] = 2.0*pnx - mom[0][i];
         mom[1][i] = 2.0*pny - mom[1][i];
         mom[2][i] = 2.0*pnz - mom[2][i];
    }
}
