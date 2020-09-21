//####### Test module for Breit-Wheeler core functions #########################

//Define Module name
 #define BOOST_TEST_MODULE "phys/breit_wheeler/core"

 #include<array>
 #include<utility>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "breit_wheeler_engine_core.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-8;
const double double_small = 1e-30;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-3;
const float float_small = 1e-20;

using namespace picsar::multi_physics::phys::breit_wheeler;

using namespace picsar::multi_physics::phys;

using namespace picsar::multi_physics::math;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

template <typename T>
T constexpr small()
{
    if(std::is_same<T,float>::value)
        return float_small;
    else
        return double_small;
}

template<typename RealType>
struct fake_T_table
{
    RealType interp(RealType chi, bool* is_out = nullptr) const {
        m_chi = chi;
        if(is_out != nullptr) *is_out = m_is_out;
        return m_res;
    }

    RealType m_res;
    bool m_is_out = false;
    mutable RealType m_chi;
};

template<typename RealType>
struct fake_P_table
{
    RealType interp(RealType chi, RealType random, bool* is_out = nullptr) const {
        m_chi = chi;
        m_random = random;
        if(is_out != nullptr) *is_out = m_is_out;
        return m_res;
    }

    RealType m_res;
    bool m_is_out = false;
    mutable RealType m_chi;
    mutable RealType m_random;
};

// ------------- Tests --------------

// ***Test Breit Wheeler optical depth

template <typename RealType>
void check_opt_depth()
{
    const auto fake_rand_zero_one_minus_epsi =
        std::array<double, 11>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999};

    for (auto ff : fake_rand_zero_one_minus_epsi){
        RealType res = get_optical_depth(ff);
        RealType expected = -log(static_cast<RealType>(1.0) - ff);
        BOOST_CHECK_EQUAL(res, expected);
    }

}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_core_opt_depth)
{
    check_opt_depth<double>();
    check_opt_depth<float>();
}

// *******************************

// ***Test Breit Wheeler dNdt

template <typename RealType, unit_system UnitSystem>
void check_dndt(RealType ref_q = one<RealType>)
{
    fake_T_table<RealType> fake_table;

    const auto chi_T_pairs = std::array<std::pair<double,double>,7>{
        std::pair<double,double>{1e-3,0.0},
        std::pair<double,double>{1e-2,3.535432626057024e-117},
        std::pair<double,double>{1e-1,5.92605753015639e-13},
        std::pair<double,double>{1,0.014135754351952334},
        std::pair<double,double>{1e1,0.10848609251601983},
        std::pair<double,double>{1e2,0.07489672111818155},
        std::pair<double,double>{1e3,0.03728739404639084}};

    const auto en_vals_SI = std::array<double,4> {
        8.187105776823886e-13,
        8.187105776823885e-12,
        8.187105776823886e-11,
        8.187105776823886e-10};

    const auto dndt_SI = std::array<std::array<double,4>,7> {
        std::array<double,4>{0.0, 0.0, 0.0, 0.0},
        std::array<double,4>{2.0029132306121232e-101, 2.002913230612124e-102, 2.0029132306121238e-103, 2.0029132306121234e-104},
        std::array<double,4>{33572.635340406545, 3357.263534040656, 335.7263534040655, 33.572635340406556},
        std::array<double,4>{8008267278956646.0, 800826727895665.0, 80082672789566.48, 8008267278956.648},
        std::array<double,4>{6.146015297711469e+17, 6.146015297711471e+16, 6146015297711470.0, 614601529771147.0},
        std::array<double,4>{4.2430912853903355e+18, 4.2430912853903366e+17, 4.243091285390336e+16, 4243091285390336.0},
        std::array<double,4>{2.112426477035044e+19, 2.1124264770350446e+18, 2.1124264770350442e+17, 2.1124264770350444e+16}};

    for (int i = 0; i < static_cast<int>(chi_T_pairs.size()); ++i){
        const auto chi = static_cast<RealType>(chi_T_pairs[i].first);
        fake_table.m_res = static_cast<RealType>(chi_T_pairs[i].second);

        for (int j = 0; j < static_cast<int>(en_vals_SI.size()); ++j ){
            const auto en = en_vals_SI[j]*conv<
                quantity::energy, unit_system::SI,
                UnitSystem, RealType>::fact(1.0, ref_q);

            const auto dndt = get_dN_dt<
                RealType,
                fake_T_table<RealType>,
                UnitSystem>(en, chi, fake_table, ref_q);

            BOOST_CHECK_EQUAL(chi,fake_table.m_chi);

            const auto sol_dndt = static_cast<RealType>(
                dndt_SI[i][j]*conv<
                    quantity::rate,unit_system::SI,
                    UnitSystem, double>::fact(1.0,ref_q));

            if(sol_dndt == static_cast<RealType>(0.0)){
                BOOST_CHECK_SMALL(dndt, small<RealType>());
            }
            else{
                BOOST_CHECK_SMALL((sol_dndt-dndt)/sol_dndt, tolerance<RealType>());
            }
        }
    }

    fake_table.m_res = 1.0;
    fake_table.m_is_out = false;
    bool is_out = false;
    get_dN_dt<RealType,fake_T_table<RealType>,UnitSystem>(
        1.0, 1.0, fake_table, ref_q, &is_out);
    BOOST_CHECK_EQUAL(is_out, false);
    fake_table.m_is_out = true;
    get_dN_dt<RealType,fake_T_table<RealType>,UnitSystem>(
        1.0, 1.0, fake_table, ref_q, &is_out);
    BOOST_CHECK_EQUAL(is_out, true);
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_core_get_dndt)
{
    const double reference_length = 800.0e-9;
    const double reference_omega = 2.0*pi<double>*light_speed<double>/
        reference_length;

    check_dndt<double, unit_system::SI>();
    check_dndt<double, unit_system::norm_omega>(reference_omega);
    check_dndt<double, unit_system::norm_lambda>(reference_length);
    check_dndt<double, unit_system::heaviside_lorentz>();
    check_dndt<float, unit_system::SI>();
    check_dndt<float, unit_system::norm_omega>(reference_omega);
    check_dndt<float, unit_system::norm_lambda>(reference_length);
    check_dndt<float, unit_system::heaviside_lorentz>();
}

// *******************************

// ***Test Breit Wheeler optical depth evolution

template <typename RealType, unit_system UnitSystem>
void check_evolve_opt_depth(RealType ref_q = one<RealType>)
{
    fake_T_table<RealType> fake_table;

    const auto chi_T_pairs = std::array<std::pair<double,double>,7>{
        std::pair<double,double>{1e-3,0.0},
        std::pair<double,double>{1e-2,3.535432626057024e-117},
        std::pair<double,double>{1e-1,5.92605753015639e-13},
        std::pair<double,double>{1,0.014135754351952334},
        std::pair<double,double>{1e1,0.10848609251601983},
        std::pair<double,double>{1e2,0.07489672111818155},
        std::pair<double,double>{1e3,0.03728739404639084}};

    const auto en_vals_SI = std::array<double,4> {
        8.187105776823886e-13,
        8.187105776823885e-12,
        8.187105776823886e-11,
        8.187105776823886e-10};

    const auto dndt_SI = std::array<std::array<double,4>,7> {
        std::array<double,4>{0.0, 0.0, 0.0, 0.0},
        std::array<double,4>{2.0029132306121232e-101, 2.002913230612124e-102, 2.0029132306121238e-103, 2.0029132306121234e-104},
        std::array<double,4>{33572.635340406545, 3357.263534040656, 335.7263534040655, 33.572635340406556},
        std::array<double,4>{8008267278956646.0, 800826727895665.0, 80082672789566.48, 8008267278956.648},
        std::array<double,4>{6.146015297711469e+17, 6.146015297711471e+16, 6146015297711470.0, 614601529771147.0},
        std::array<double,4>{4.2430912853903355e+18, 4.2430912853903366e+17, 4.243091285390336e+16, 4243091285390336.0},
        std::array<double,4>{2.112426477035044e+19, 2.1124264770350446e+18, 2.1124264770350442e+17, 2.1124264770350444e+16}};

    const auto dt_SI =  std::array<double,4>{1e-18, 1e-15, 1e-12, 1e-9};

    for (int i = 0; i < static_cast<int>(chi_T_pairs.size()); ++i){
        const auto chi = static_cast<RealType>(chi_T_pairs[i].first);
        fake_table.m_res = static_cast<RealType>(chi_T_pairs[i].second);

        for (int j = 0; j < static_cast<int>(en_vals_SI.size()); ++j ){
            const auto en = en_vals_SI[j]*conv<
                quantity::energy, unit_system::SI,
                UnitSystem, RealType>::fact(1.0, ref_q);

            for (auto t_dt : dt_SI){
                const RealType init_opt = 1.0;
                auto opt_depth = init_opt;

                const auto dt = t_dt*conv<
                    quantity::time, unit_system::SI,
                    UnitSystem, RealType>::fact(1.0, ref_q);

                fake_table.m_is_out = false;
                const bool in_table_flag = evolve_optical_depth<
                    RealType,
                    fake_T_table<RealType>,
                    UnitSystem>(en, chi, dt, opt_depth, fake_table, ref_q);

                BOOST_CHECK_EQUAL(chi,fake_table.m_chi);
                BOOST_CHECK_EQUAL(in_table_flag, true);

                const auto sol_dndt = static_cast<RealType>(
                    dndt_SI[i][j]*conv<
                        quantity::rate,unit_system::SI,
                        UnitSystem, double>::fact(1.0,ref_q));

                const RealType sol = init_opt - sol_dndt*dt;

                BOOST_CHECK_SMALL((sol-opt_depth)/sol, tolerance<RealType>());

                fake_table.m_is_out = true;
                const bool in_table_flag_2 = evolve_optical_depth<
                    RealType,
                    fake_T_table<RealType>,
                    UnitSystem>(en, chi, dt, opt_depth, fake_table, ref_q);
                BOOST_CHECK_EQUAL(in_table_flag_2, false);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_core_evolve_opt_depth)
{
    const double reference_length = 800.0e-9;
    const double reference_omega = 2.0*pi<double>*light_speed<double>/
        reference_length;

    check_evolve_opt_depth<double, unit_system::SI>();
    check_evolve_opt_depth<double, unit_system::norm_omega>(reference_omega);
    check_evolve_opt_depth<double, unit_system::norm_lambda>(reference_length);
    check_evolve_opt_depth<double, unit_system::heaviside_lorentz>();
    check_evolve_opt_depth<float, unit_system::SI>();
    check_evolve_opt_depth<float, unit_system::norm_omega>(reference_omega);
    check_evolve_opt_depth<float, unit_system::norm_lambda>(reference_length);
    check_evolve_opt_depth<float, unit_system::heaviside_lorentz>();
}

// *******************************

// ***Test Breit Wheeler pair production

template <typename RealType>
void aux_check_small_or_zero(vec3<RealType> res, vec3<RealType> expected)
{
    for(int i = 0; i < 3; ++i){
        if(expected[i] == zero<RealType>){
            BOOST_CHECK_SMALL(expected[i]-res[i], small<RealType>());
        }
        else{
            BOOST_CHECK_SMALL((expected[i]-res[i])/res[i], tolerance<RealType>());
        }
    }
}

template <typename RealType, unit_system UnitSystem>
void check_pair_production(RealType ref_q = one<RealType>)
{
    fake_P_table<RealType> fake_table;

    const auto chi_photons = std::array<double,7>{
        1.0e-3,1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3};

    const auto dirs = std::array<vec3<double>,4>{
        vec3<double>{1.0, 0.0, 0.0},
        vec3<double>{0.0, 1.0, 0.0},
        vec3<double>{0.0, 0.0, 1.0},
        vec3<double>{1.0, 1.0, 1.0}/sqrt(3.0)};

    const auto photon_momenta_SI = std::array<double,3> {
        2.7309245307378235e-21,
        2.730924530737823e-20,
        2.730924530737823e-19};

    const auto random = std::array<double,4> {
        0.0, 0.3,0.9,0.999};

    const auto chi_ele_fracs = std::array<double, 7> {
        1.0e-2, 1.0e-1, 0.2, 0.5, 0.7, 0.9, 0.99};

    for (auto tchi_phot : chi_photons){
        for (auto dir : dirs){
            for(auto mom : photon_momenta_SI){
                for (auto rand : random){
                    for (auto chi_ele_frac : chi_ele_fracs){
                        const RealType chi_photon = tchi_phot;
                        const RealType chi_ele = chi_ele_frac*tchi_phot;
                        fake_table.m_res = chi_ele;
                        fake_table.m_is_out = false;

                        const RealType pmom = static_cast<RealType>(mom);

                        const auto photon_momentum =
                            vec3<RealType>{
                                static_cast<RealType>(dir[0]),
                                static_cast<RealType>(dir[1]),
                                static_cast<RealType>(dir[2])} *
                                pmom *
                            conv<quantity::momentum,
                                unit_system::SI,UnitSystem,
                                RealType>::fact(1.0, ref_q);

                        const auto rand_num = static_cast<RealType>(rand);

                        vec3<RealType> mom_ele;
                        vec3<RealType> mom_pos;

                        const bool is_in_table = generate_breit_wheeler_pairs<
                            RealType,
                            fake_P_table<RealType>,
                            UnitSystem>(
                                chi_photon, photon_momentum, rand_num,
                                fake_table, mom_ele, mom_pos, ref_q);

                        BOOST_CHECK_EQUAL(fake_table.m_random,rand_num);
                        BOOST_CHECK_EQUAL(fake_table.m_chi,chi_photon);
                        BOOST_CHECK_EQUAL(!fake_table.m_is_out, is_in_table);

                        const auto me_c = electron_mass<>*light_speed<>;
                        const auto gamma_phot = mom/(me_c);
                        const auto en_ele = (1. + (gamma_phot-2.)*chi_ele/tchi_phot);
                        const auto en_pos = (1. + (gamma_phot-2.)*(tchi_phot - chi_ele)/tchi_phot);
                        const auto t_mom_ele = me_c*sqrt( en_ele*en_ele - 1);
                        const auto t_mom_pos = me_c*sqrt( en_pos*en_pos - 1);

                        const auto tv_mom_ele = t_mom_ele * dir *
                            conv<quantity::momentum,
                                unit_system::SI,UnitSystem,
                                double>::fact(1.0, ref_q);
                        const auto tv_mom_pos = t_mom_pos * dir *
                            conv<quantity::momentum,
                                unit_system::SI,UnitSystem,
                                double>::fact(1.0, ref_q);

                        const auto exp_mom_ele = vec3<RealType>{
                            static_cast<RealType>(tv_mom_ele[0]),
                            static_cast<RealType>(tv_mom_ele[1]),
                            static_cast<RealType>(tv_mom_ele[2])};

                        const auto exp_mom_pos = vec3<RealType>{
                            static_cast<RealType>(tv_mom_pos[0]),
                            static_cast<RealType>(tv_mom_pos[1]),
                            static_cast<RealType>(tv_mom_pos[2])};

                        aux_check_small_or_zero(mom_ele, exp_mom_ele);
                        aux_check_small_or_zero(mom_pos, exp_mom_pos);

                        fake_table.m_is_out = true;
                        const bool is_in_table_2 = generate_breit_wheeler_pairs<
                            RealType,
                            fake_P_table<RealType>,
                            UnitSystem>(
                                chi_photon, photon_momentum, rand_num,
                                fake_table, mom_ele, mom_pos, ref_q);
                        BOOST_CHECK_EQUAL(!fake_table.m_is_out, is_in_table_2);
                    }
                }
            }
        }


    }
}

BOOST_AUTO_TEST_CASE( picsar_breit_wheeler_core_pair_production)
{
    const double reference_length = 800.0e-9;
    const double reference_omega = 2.0*pi<double>*light_speed<double>/
        reference_length;

    check_pair_production<double, unit_system::SI>();
    check_pair_production<double, unit_system::norm_omega>(reference_omega);
    check_pair_production<double, unit_system::norm_lambda>(reference_length);
    check_pair_production<double, unit_system::heaviside_lorentz>();
    check_pair_production<float, unit_system::SI>();
    check_pair_production<float, unit_system::norm_omega>(reference_omega);
    check_pair_production<float, unit_system::norm_lambda>(reference_length);
    check_pair_production<float, unit_system::heaviside_lorentz>();
}

// *******************************
