#include "landau_lifshitz.h"

using namespace std;
using namespace picsar::multi_physics;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING CHECK MASS NORMALIZATION CONVENTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void  picsar::multi_physics::boris_plus_landau_lifshitz_push(
vector<double>& px, vector<double>& py, vector<double>& pz,
const vector<double>& ex, const vector<double>& ey, const vector<double>& ez,
const vector<double>& bx, const vector<double>& by, const vector<double>& bz,
double mass, double charge, double dt, double lambda){

    if(mass == 0 || charge == 0 || lambda == 0 || dt == 0)
        return;

    double beta = 0.5 * charge * dt / mass;
    double inv_dt = 1.0/dt;

//Radiation reaction coefficient
    double rr_coeff = (4.0 * M_PI * classical_electron_radius / 3.0 / lambda);

    double old_px, old_py, old_pz;

    int num_particles = px.size();
    for (int i = 0; i < num_particles; i++){
        old_px = px[i];
        old_py = py[i];
        old_pz = pz[i];

        double ptx = px[i] + beta *  ex[i];
        double pty = py[i] + beta *  ey[i];
        double ptz = pz[i] + beta *  ez[i];

        double inv_gamman_times_beta = beta/sqrt(1.0 + (px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i])/(mass*mass));

        double bbx = inv_gamman_times_beta * bx[i];
        double bby = inv_gamman_times_beta * by[i];
        double bbz = inv_gamman_times_beta * bz[i];

        double bb2 = bbx*bbx + bby*bby + bbz*bbz;
        double coeff = 1.0 / (1.0 + bb2);

        double ptvbx = pty*bbz - ptz*bby;
        double ptvby = ptz*bbx - ptx*bbz;
        double ptvbz = ptx*bby - pty*bbx;

        double ptdb = ptx * bbx + pty * bby + ptz * bbz;

        double bptdbx = ptdb*bbx;
        double bptdby = ptdb*bby;
        double bptdbz = ptdb*bbz;

        double pex = coeff*(ptx + ptvbx + bptdbx);
        double pey = coeff*(pty + ptvby + bptdby);
        double pez = coeff*(ptz + ptvbz + bptdbz);

        px[i] = 2.0*pex - px[i];
        py[i] = 2.0*pey - py[i];
        pz[i] = 2.0*pez - pz[i];

        //Momenta at time n
        double pxn = (px[i] + old_px)*0.5;
        double pyn = (py[i] + old_py)*0.5;
        double pzn = (pz[i] + old_pz)*0.5;

        //Inverse gamma at time n
        double gn =  sqrt(1.0 + (pxn*pxn + pyn*pyn + pzn*pzn)/(mass*mass));
        double inv_gn = 1.0/gn;

        //Velocities at time n
        double vxn = pxn*inv_gn;
        double vyn = pyn*inv_gn;
        double vzn = pzn*inv_gn;

        //Lorentz force
        double flx = (px[i] - old_px)*inv_dt;
        double fly = (py[i] - old_py)*inv_dt;
        double flz = (pz[i] - old_pz)*inv_dt;

        //RR force (M.Tamburini scheme)
        double fl2 = flx*flx + fly*fly + flz*flz;
        double vdote = vxn*ex[i] +  vyn*ey[i] +  vzn*ez[i];
        double vdote2 = vdote*vdote;

        double mul = -rr_coeff*(gn*gn)*(fl2 - vdote2)*dt;
        double prrx = mul*vxn;
        double prry = mul*vyn;
        double prrz = mul*vzn;

        px[i] += prrx;
        py[i] += prry;
        pz[i] += prrz;

    }
}
