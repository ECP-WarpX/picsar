#include "chi_calc_functions.h"

using namespace std;
using namespace picsar::multi_physics;

void picsar::multi_physics::chi_lepton_lambda_vec(
    vector<double>& px, vector<double>& py, vector<double>& pz,
    const vector<double>& ex, const vector<double>& ey, const vector<double>& ez,
    const vector<double>& bx, const vector<double>& by, const vector<double>& bz,
    double mass, double lambda,
    vector<double>& chi){

        int num_part = px.size();
        for(int i = 0; i < num_part; i++){
            chi[i] = chi_lepton_lambda({px[i], py[i], pz[i]},
                {ex[i], ey[i], ez[i], bx[i], by[i], bz[i]},
                mass, lambda);
        }

}

void picsar::multi_physics::chi_photon_lambda_vec(
    vector<double>& px, vector<double>& py, vector<double>& pz,
    const vector<double>& ex, const vector<double>& ey, const vector<double>& ez,
    const vector<double>& bx, const vector<double>& by, const vector<double>& bz,
    double lambda,
    vector<double>& chi){
        int num_part = px.size();

        for(int i = 0; i < num_part; i++){
            chi[i] = chi_photon_lambda({px[i], py[i], pz[i]},
                {ex[i], ey[i], ez[i], bx[i], by[i], bz[i]},
                lambda);
        }

}
