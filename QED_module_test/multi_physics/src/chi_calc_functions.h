#ifndef __PMP_CHICALC__
#define __PMP_CHICALC__

#include <array>
#include <vector>
#include <cmath>

#include "commons.h"

namespace picsar{
    namespace multi_physics{

        inline double chi_lepton_inv_schwinger(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double mass, double normalized_inv_schwinger_field);
        inline double chi_lepton_lambda(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double mass,  double lambda);
        inline double chi_photon_inv_schwinger(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double normalized_inv_schwinger_field);
        inline double chi_photon_lambda(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double lambda);

        void chi_lepton_lambda_vec(
            std::vector<double>& px, std::vector<double>& py, std::vector<double>& pz,
            const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez,
            const std::vector<double>& bx, const std::vector<double>& by, const std::vector<double>& bz,
            double mass, double lambda,
            std::vector<double>& chi
        );

        void chi_photon_lambda_vec(
            std::vector<double>& px, std::vector<double>& py, std::vector<double>& pz,
            const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez,
            const std::vector<double>& bx, const std::vector<double>& by, const std::vector<double>& bz,
            double lambda,
            std::vector<double>& chi
        );

        inline double chi_lepton_inv_schwinger(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double mass, double normalized_inv_schwinger_field){
            double gamma = sqrt(1.0 + mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
            double inv_gamma = 1.0/gamma;

            double vx = inv_gamma * mom[0];
            double vy = inv_gamma * mom[1];
            double vz = inv_gamma * mom[2];

            double vdotex =  vx * emfl[0];
            double vdotey =  vy * emfl[1];
            double vdotez =  vz * emfl[2];

            double vdote2 = vdotex*vdotex + vdotey*vdotey + vdotez*vdotez;

            double eplusvcrossbx = emfl[0] + (vy * emfl[5] - vz*emfl[4]);
            double eplusvcrossby = emfl[1] + (vz * emfl[3] - vx*emfl[5]);
            double eplusvcrossbz = emfl[2] + (vx * emfl[4] - vy*emfl[3]);

            double eplusvcrossb2 = eplusvcrossbx*eplusvcrossbx + eplusvcrossby*eplusvcrossby + eplusvcrossbz*eplusvcrossbz;

            return gamma * mass * sqrt(fabs(vdote2 - eplusvcrossb2))*normalized_inv_schwinger_field;

        }

        inline double chi_lepton_lambda(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double mass, double lambda){
            return chi_lepton_inv_schwinger(mom, emfl, 1.0/picsar::multi_physics::calc_schwinger_given_lambda(lambda), mass);
        }

        inline double chi_photon_inv_schwinger(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double normalized_inv_schwinger_field){
            double p = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
            double inv_p = 1.0/p;

            double cx = mom[0]*inv_p;
            double cy = mom[1]*inv_p;
            double cz = mom[2]*inv_p;

            double cdote = cx*emfl[0] + cy*emfl[1] + cz*emfl[2];
            double eperx = cdote*cx;
            double epery = cdote*cy;
            double eperz = cdote*cz;

            double ccrossbx = cy * emfl[5] - cz*emfl[4];
            double ccrossby = cz * emfl[3] - cx*emfl[5];
            double ccrossbz = cx * emfl[4] - cy*emfl[3];

            double eplusccrossbx = eperx + ccrossbx;
            double eplusccrossby = epery + ccrossby;
            double eplusccrossbz = eperz + ccrossbz;

            double eplusccrossbz2 = eplusccrossbx*eplusccrossbx + eplusccrossby*eplusccrossby + eplusccrossbz*eplusccrossbz;

            return p * sqrt(eplusccrossbz2)*normalized_inv_schwinger_field;
        }

        inline double chi_photon_lambda(const std::array<double, 3>& mom, const std::array<double, 6>& emfl, double lambda){
            return chi_photon_inv_schwinger(mom, emfl, 1.0/picsar::multi_physics::calc_schwinger_given_lambda(lambda));
        }

    }
}

#endif
