#include "leptons.h"

using namespace testbed;

leptons::leptons(std::string name):species(name){}

leptons::~leptons(){}

void leptons::push_positions(ttime dt){
    int num_particles = pos[0].size();

    for (int i = 0; i < num_particles; i++){
        double inv_ptotal = 1.0/sqrt(1.0 + mom[0][i]*mom[0][i] + mom[1][i]*mom[1][i] + mom[2][i]*mom[2][i]);
        double coeff = dt*inv_ptotal;

        pos[0][i] += coeff*mom[0][i];
        pos[1][i] += coeff*mom[1][i];
        pos[2][i] += coeff*mom[2][i];
    }
}

void leptons::push_momenta(ttime dt){
    if (!is_boris_replaced)
        do_boris_pusher(dt);
    else
        do_alternative_push_momenta(dt);
}

//Boris pusher
void leptons::do_boris_pusher(double dt){
    double beta = 0.5 * charge * dt / mass;

    int num_particles = pos[0].size();
    for (int i = 0; i < num_particles; i++){

        double ptx = mom[0][i] + beta *  fields[0][i];
        double pty = mom[1][i] + beta *  fields[1][i];
        double ptz = mom[2][i] + beta *  fields[2][i];

        double inv_gamman_times_beta = beta/sqrt(1.0 + mom[0][i]*mom[0][i] + mom[1][i]*mom[1][i] + mom[2][i]*mom[2][i]);

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

        double pnx = coeff*(ptx + ptvbx + bptdbx);
        double pny = coeff*(pty + ptvby + bptdby);
        double pnz = coeff*(ptz + ptvbz + bptdbz);

         mom[0][i] = 2.0*pnx - mom[0][i];
         mom[1][i] = 2.0*pny - mom[1][i];
         mom[2][i] = 2.0*pnz - mom[2][i];
    }
}

void leptons::replace_pusher_momenta(mom_pusher_function pusher){
    is_boris_replaced = true;
    this->pusher = pusher;
}

void leptons::do_alternative_push_momenta(ttime dt){
        pusher(mom, fields, mass, charge, dt);
}
