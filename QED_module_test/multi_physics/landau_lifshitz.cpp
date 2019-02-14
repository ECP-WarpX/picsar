#include "landau_lifshitz.h"

using namespace picsar;

void picsar::boris_plus_landau_lifshitz_push(momenta_list& mom, const em_field_list& fields, double mass, double charge, double dt, double lambda){
    double beta = 0.5 * charge * dt / mass;
    double inv_dt = 1.0/dt;

//Radiation reaction coefficient
    double rr_coeff = (4.0 * M_PI * classical_electron_radius / 3.0 / lambda);

    double old_px, old_py, old_pz;

    int num_particles = mom[0].size();
    for (int i = 0; i < num_particles; i++){
        old_px = mom[0][i];
        old_py = mom[1][i];
        old_pz = mom[2][i];

        double ptx = mom[0][i] + beta *  fields[0][i];
        double pty = mom[1][i] + beta *  fields[1][i];
        double ptz = mom[2][i] + beta *  fields[2][i];

        double gamman = sqrt(mass*mass + mom[0][i]*mom[0][i] + mom[1][i]*mom[1][i] + mom[2][i]*mom[2][i]);
        double inv_gamman_times_beta = beta/gamman;

        double bx = inv_gamman_times_beta * fields[3][i];
        double by = inv_gamman_times_beta * fields[4][i];
        double bz = inv_gamman_times_beta * fields[5][i];

        double b2 = bx*bx + by*by + bz*bz;
        double coeff = 1.0 / (1.0 + b2);

        double ptvbx = pty*bz - ptz*by;
        double ptvby = ptz*bx - ptx*bz;
        double ptvbz = ptx*by - pty*bx;

        double ptdb = ptx * bx + pty * by + ptz * bz;

        double bptdbx = ptdb*bx;
        double bptdby = ptdb*by;
        double bptdbz = ptdb*bz;

        double pex = coeff*(ptx + ptvbx + bptdbx);
        double pey = coeff*(pty + ptvby + bptdby);
        double pez = coeff*(ptz + ptvbz + bptdbz);

        mom[0][i] = 2.0*pex - mom[0][i];
        mom[1][i] = 2.0*pey - mom[1][i];
        mom[2][i] = 2.0*pez - mom[2][i];

        //Momenta at time n
        double pxn = (mom[0][i] + old_px)*0.5;
        double pyn = (mom[1][i] + old_py)*0.5;
        double pzn = (mom[2][i] + old_pz)*0.5;

        //Inverse gamma at time n
        double inv_gn = sqrt(mass*mass + pxn*pxn + pyn*pyn + pzn*pzn);

        //Velocities at time n
        double vxn = pxn*inv_gn;
        double vyn = pyn*inv_gn;
        double vzn = pzn*inv_gn;

        //Lorentz force
        double flx = (mom[0][i] - old_px)*inv_dt;
        double fly = (mom[1][i] - old_py)*inv_dt;
        double flz = (mom[2][i] - old_pz)*inv_dt;

        //RR force (M.Tamburini scheme)
        double fl2 = flx*flx + fly*fly + flz*flz;
        double vdote = vxn*fields[0][i] +  vyn*fields[1][i] +  vzn*fields[2][i];
        double vdote2 = vdote*vdote;

        double mul = -rr_coeff*(gamman*gamman)*rr_coeff*(fl2 - vdote2)*dt;
        double prrx = mul*vxn;
        double prry = mul*vyn;
        double prrz = mul*vzn;

        mom[0][i] += prrx;
        mom[1][i] += prry;
        mom[2][i] += prrz;

    }
}
