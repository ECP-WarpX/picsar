//####### Test module for Quantum Synchrotron core functions ###################

//Define Module name
 #define BOOST_TEST_MODULE "phys/quantum_sync/core"

 #include<array>
 #include<utility>

//Include Boost unit tests library & library for floating point comparison
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include "quantum_sync_engine_core.hpp"

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-8;
const double double_small = 1e-30;

//Tolerance for single precision calculations
const float float_tolerance = 3.0e-3;
const float float_small = 1e-20;

using namespace picsar::multi_physics::phys::quantum_sync;

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
struct fake_G_table
{
    RealType interp(RealType chi, bool* is_out = nullptr) const {
        m_chi = chi;
        if(is_out != nullptr) *is_out = m_is_out;
        return static_cast<RealType>(m_res);
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
        return static_cast<RealType>(m_res);
    }

    RealType m_res;
    bool m_is_out = false;
    mutable RealType m_chi;
    mutable RealType m_random;
};

// ------------- Tests --------------

// ***Test Quantum Synchrotron optical depth

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

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_core_opt_depth)
{
    check_opt_depth<double>();
    check_opt_depth<float>();
}

// *******************************

// ***Test Quantum Synchrotron dNdt

template <typename RealType, unit_system UnitSystem>
void check_dndt(RealType ref_q = one<RealType>)
{
    fake_G_table<RealType> fake_table;

    const auto chi_G_pairs = std::array<std::pair<double,double>,7>{
            std::pair<double,double>{1e-3,0.0021630710451102635},
            std::pair<double,double>{1e-2,0.02145778613250966},
            std::pair<double,double>{1e-1,0.20141650057288696},
            std::pair<double,double>{1,1.5508709239783094},
            std::pair<double,double>{1e1,9.170292626506058},
            std::pair<double,double>{1e2,46.02341774244706},
            std::pair<double,double>{1e3,217.8438638968691}};

    const auto en_vals_SI = std::array<double,4> {
        8.187105776823886e-13,
        8.187105776823885e-12,
        8.187105776823886e-11,
        8.187105776823886e-10};

    const auto dndt_SI = std::array<std::array<double,4>,7> {
        std::array<double,4>{816956805243417.9, 81695680524341.8,
            8169568052434.179, 816956805243.4178},
        std::array<double,4>{8104257345610210.0, 810425734561021.0,
            81042573456102.1, 8104257345610.21},
        std::array<double,4>{7.607174124183555e+16, 7607174124183556.0,
            760717412418355.6, 76071741241835.56},
        std::array<double,4>{5.857387616843818e+17, 5.857387616843818e+16,
            5857387616843818.0, 585738761684381.8},
        std::array<double,4>{3.463470598542345e+18, 3.4634705985423456e+17,
            3.4634705985423456e+16, 3463470598542345.5},
        std::array<double,4>{1.7382297456316862e+19, 1.7382297456316861e+18,
            1.738229745631686e+17, 1.7382297456316862e+16},
        std::array<double,4>{8.227608958724519e+19, 8.227608958724519e+18,
            8.227608958724518e+17, 8.227608958724518e+16}};

    for (int i = 0; i < chi_G_pairs.size(); ++i){
        const auto chi = static_cast<RealType>(chi_G_pairs[i].first);
        fake_table.m_res = static_cast<RealType>(chi_G_pairs[i].second);

        for (int j = 0; j < en_vals_SI.size(); ++j ){
            const auto en = en_vals_SI[j]*conv<
                quantity::energy, unit_system::SI,
                UnitSystem, RealType>::fact(1.0, ref_q);

            const auto dndt = get_dN_dt<
                RealType,
                fake_G_table<RealType>,
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
    get_dN_dt<RealType,fake_G_table<RealType>,UnitSystem>(
        1.0, 1.0, fake_table, ref_q, &is_out);
    BOOST_CHECK_EQUAL(is_out, false);
    fake_table.m_is_out = true;
    get_dN_dt<RealType,fake_G_table<RealType>,UnitSystem>(
        1.0, 1.0, fake_table, ref_q, &is_out);
    BOOST_CHECK_EQUAL(is_out, true);
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_core_get_dndt)
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

// ***Test Quantum Synchrotron optical depth evolution

template <typename RealType, unit_system UnitSystem>
void check_evolve_opt_depth(RealType ref_q = one<RealType>)
{
    fake_G_table<RealType> fake_table;

    const auto chi_G_pairs = std::array<std::pair<double,double>,7>{
        std::pair<double,double>{1e-3,0.0021630710451102635},
        std::pair<double,double>{1e-2,0.02145778613250966},
        std::pair<double,double>{1e-1,0.20141650057288696},
        std::pair<double,double>{1,1.5508709239783094},
        std::pair<double,double>{1e1,9.170292626506058},
        std::pair<double,double>{1e2,46.02341774244706},
        std::pair<double,double>{1e3,217.8438638968691}};

    const auto en_vals_SI = std::array<double,4> {
        8.187105776823886e-13,
        8.187105776823885e-12,
        8.187105776823886e-11,
        8.187105776823886e-10};

    const auto dndt_SI = std::array<std::array<double,4>,7> {
        std::array<double,4>{816956805243417.9, 81695680524341.8, 8169568052434.179, 816956805243.4178},
        std::array<double,4>{8104257345610210.0, 810425734561021.1, 81042573456102.1, 8104257345610.21},
        std::array<double,4>{7.607174124183555e+16, 7607174124183557.0, 760717412418355.6, 76071741241835.56},
        std::array<double,4>{5.857387616843818e+17, 5.857387616843819e+16, 5857387616843818.0, 585738761684381.8},
        std::array<double,4>{3.463470598542345e+18, 3.4634705985423456e+17, 3.4634705985423456e+16, 3463470598542345.5},
        std::array<double,4>{1.7382297456316862e+19, 1.7382297456316864e+18, 1.738229745631686e+17, 1.7382297456316862e+16},
        std::array<double,4>{8.227608958724519e+19, 8.22760895872452e+18,8.227608958724518e+17, 8.227608958724518e+16}};

    const auto dt_SI =  std::array<double,4>{1e-18, 1e-15, 1e-12, 1e-9};

    for (int i = 0; i < chi_G_pairs.size(); ++i){
        const auto chi = static_cast<RealType>(chi_G_pairs[i].first);
        fake_table.m_res = static_cast<RealType>(chi_G_pairs[i].second);

        for (int j = 0; j < en_vals_SI.size(); ++j ){
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
                    fake_G_table<RealType>,
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
                    fake_G_table<RealType>,
                    UnitSystem>(en, chi, dt, opt_depth, fake_table, ref_q);

                BOOST_CHECK_EQUAL(in_table_flag_2, false);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_core_evolve_opt_depth)
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

// ***Test Quantum Synchrotron photon emission

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
void check_photon_emission(RealType ref_q = one<RealType>)
{
    fake_P_table<RealType> fake_table;

    const auto chi_eles = std::array<double,7>{
        1.0e-3,1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3};

    const auto dirs = std::array<vec3<double>,4>{
        vec3<double>{1.0, 0.0, 0.0},
        vec3<double>{0.0, 1.0, 0.0},
        vec3<double>{0.0, 0.0, 1.0},
        vec3<double>{1.0, 1.0, 1.0}/sqrt(3.0)};

    const auto ele_momenta_SI = std::array<double,3> {
        2.7309245307378235e-21,
        2.730924530737823e-20,
        2.730924530737823e-19};

    const auto random = std::array<double,4> {
        0.0, 0.3,0.9,0.999};

    const auto chi_phot_fracs = std::array<double, 7> {
        1.0e-2, 1.0e-1, 0.2, 0.5, 0.7, 0.9, 0.99};

    for (auto t_chi_ele : chi_eles){
        for (auto dir : dirs){
            for(auto tmom : ele_momenta_SI){
                for (auto rand : random){
                    for (auto chi_phot_frac : chi_phot_fracs){
                        const auto chi_ele = static_cast<RealType>(t_chi_ele);
                        const RealType chi_phot = chi_phot_frac*chi_ele;
                        fake_table.m_res = chi_phot;
                        fake_table.m_is_out = false;

                        const auto mom = static_cast<RealType>(tmom);

                        auto ele_momentum =
                            vec3<RealType>{
                                static_cast<RealType>(dir[0]),
                                static_cast<RealType>(dir[1]),
                                static_cast<RealType>(dir[2])} *
                                mom *
                            conv<quantity::momentum,
                                unit_system::SI,UnitSystem,
                                RealType>::fact(1.0, ref_q);

                        const auto rand_num = static_cast<RealType>(rand);

                        vec3<RealType> mom_phot;
                        const bool is_in_table = generate_photon_update_momentum<
                            RealType,
                            fake_P_table<RealType>,
                            UnitSystem>(
                                chi_ele, ele_momentum, rand_num,
                                fake_table, mom_phot, ref_q);

                        BOOST_CHECK_EQUAL(fake_table.m_random, rand_num);
                        BOOST_CHECK_EQUAL(fake_table.m_chi, chi_ele);
                        BOOST_CHECK_EQUAL(!fake_table.m_is_out, is_in_table);

                        const auto me_c = electron_mass<>*light_speed<>;
                        const auto gamma_ele = sqrt(1.0 + (mom/me_c)*(mom/me_c));
                        const auto gamma_phot = (gamma_ele-1.)*(chi_phot/chi_ele);
                        const auto t_mom_phot = me_c*gamma_phot;

                        const auto tv_mom_phot = t_mom_phot * dir *
                            conv<quantity::momentum,
                                unit_system::SI,UnitSystem,
                                double>::fact(1.0, ref_q);

                        const auto exp_mom_phot = vec3<RealType>{
                            static_cast<RealType>(tv_mom_phot[0]),
                            static_cast<RealType>(tv_mom_phot[1]),
                            static_cast<RealType>(tv_mom_phot[2])};

                        aux_check_small_or_zero(mom_phot, exp_mom_phot);

                        fake_table.m_is_out = true;
                        const bool is_in_table_2 = generate_photon_update_momentum<
                            RealType,
                            fake_P_table<RealType>,
                            UnitSystem>(
                                chi_ele, ele_momentum, rand_num,
                                fake_table, mom_phot, ref_q);
                        BOOST_CHECK_EQUAL(!fake_table.m_is_out, is_in_table_2);
                    }
                }
            }
        }


    }
}

BOOST_AUTO_TEST_CASE( picsar_quantum_sync_photon_emission)
{
    const double reference_length = 800.0e-9;
    const double reference_omega = 2.0*pi<double>*light_speed<double>/
        reference_length;

    check_photon_emission<double, unit_system::SI>();
    check_photon_emission<double, unit_system::norm_omega>(reference_omega);
    check_photon_emission<double, unit_system::norm_lambda>(reference_length);
    check_photon_emission<double, unit_system::heaviside_lorentz>();
    check_photon_emission<float, unit_system::SI>();
    check_photon_emission<float, unit_system::norm_omega>(reference_omega);
    check_photon_emission<float, unit_system::norm_lambda>(reference_length);
    check_photon_emission<float, unit_system::heaviside_lorentz>();
}

// *******************************
