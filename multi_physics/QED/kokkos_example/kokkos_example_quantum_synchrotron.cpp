#include "kokkos_example_commons.hpp"

//Parameters of the test case
const unsigned int how_many_particles = 10'000'000;
const unsigned int how_many_repetitions = 1;
const double dt_test= 1e-18;
const double table_chi_min = 0.01;
const double table_frac_min = 1.0e-12;
const double table_chi_max = 1000.0;
const int table_chi_size = 128;
const int table_frac_size = 128;
const int random_seed = 22051988;
const double E_min = -0.01*Es;
const double E_max = 0.01*Es;
const double B_min = E_min/pxr::light_speed<>;
const double B_max = E_max/pxr::light_speed<>;
const double P_min = -100*mec<>;
const double P_max = 100*mec<>;
//__________________________________________________

template <typename Real, typename Vector>
auto generate_dndt_table(Real chi_min, Real chi_max, int chi_size)
{
    std::cout << "Preparing dndt table [" << typeid(Real).name() << ", " << chi_size <<"]...\n";
    std::cout.flush();

    pxr_qs::dndt_lookup_table_params<Real> qs_params{chi_min, chi_max, chi_size};

	auto table = pxr_qs::dndt_lookup_table<
        Real, Vector>{qs_params};

    table.generate();

    return table;
}

template <typename Real, typename Vector>
auto generate_photon_emission_table(
    Real chi_min, Real chi_max, Real frac_min, int chi_size, int frac_size)
{
    std::cout << "Preparing photon emission table [" << typeid(Real).name() << ", " << chi_size << " x " << frac_size <<"]...\n";
    std::cout.flush();

    pxr_qs::photon_emission_lookup_table_params<Real> qs_params{
        chi_min, chi_max, frac_min, chi_size, frac_size};

	auto table = pxr_qs::photon_emission_lookup_table<
        Real, Vector>{qs_params};

    table.template generate();

    return table;
}

template <typename Real>
bool fill_opt_test(
    ParticleData<Real>& pdata,
    const int repetitions,
    Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    for(int rr = 0; rr < repetitions; ++rr){
        const auto num_particles = pdata.num_particles;
        Kokkos::parallel_for("FillOpt_"+get_type_name<Real>(),
            num_particles, KOKKOS_LAMBDA(int i){
            auto rand_gen = rand_pool.get_state();
            pdata.m_fields.opt(i) =
                pxr_qs::get_optical_depth<Real>(
                    get_rand<Real>::get(rand_gen));
            rand_pool.free_state(rand_gen);
        });
    }

    return check(pdata.m_fields.opt, true, false);
}


template <typename Real, typename TableType>
bool evolve_optical_depth(
    ParticleData<Real>& pdata, const TableType& ref_table, Real dt, const int repetitions)
{
    for(int rr = 0; rr < repetitions; ++rr){
        const auto num_particles = pdata.num_particles;
        Kokkos::parallel_for("EvolveOpt", num_particles, KOKKOS_LAMBDA(int i){
            const auto px = pdata.m_momentum(i,0);
            const auto py = pdata.m_momentum(i,1);
            const auto pz = pdata.m_momentum(i,2);
            const auto ex = pdata.m_fields.Ex(i);
            const auto ey = pdata.m_fields.Ey(i);
            const auto ez = pdata.m_fields.Ez(i);
            const auto bx = pdata.m_fields.Bx(i);
            const auto by = pdata.m_fields.By(i);
            const auto bz = pdata.m_fields.Bz(i);
            auto& opt = pdata.m_fields.opt(i);
            const auto ee =
                pxr::compute_gamma_ele_pos<Real>(px, py, pz)*mec2<Real>;
            const auto chi =
                pxr::chi_ele_pos<Real, pxr::unit_system::SI>(
                    px, py ,pz, ex, ey, ez, bx, by, bz);
            pxr_qs::evolve_optical_depth<Real, TableType>(
                ee, chi, dt, opt, ref_table);
        });
    }

    return check(pdata.m_fields.opt, true, false);
}


template <typename Real, typename TableType>
bool generate_photons(
    ParticleData<Real>& pdata, const TableType& ref_table, const int repetitions,
    Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    const auto num_particles = pdata.num_particles;

    auto photon_momentum = init_multi_comp_view_with_random_content<Real>(
        "photon_momentum", 0.0, 0.0, num_particles, rand_pool);

    for(int rr = 0; rr < repetitions; ++rr){
        Kokkos::parallel_for("PhotEm", num_particles, KOKKOS_LAMBDA(int i){
            const auto px = pdata.m_momentum(i,0);
            const auto py = pdata.m_momentum(i,1);
            const auto pz = pdata.m_momentum(i,2);
            const auto ex = pdata.m_fields.Ex(i);
            const auto ey = pdata.m_fields.Ey(i);
            const auto ez = pdata.m_fields.Ez(i);
            const auto bx = pdata.m_fields.Bx(i);
            const auto by = pdata.m_fields.By(i);
            const auto bz = pdata.m_fields.Bz(i);

            const auto chi = pxr::chi_ele_pos<Real, pxr::unit_system::SI>(
                px, py ,pz, ex, ey, ez, bx, by, bz);

            auto rand_gen = rand_pool.get_state();

            auto e_phot = pxr_m::vec3<Real>{0,0,0};

            auto p_mom = pxr_m::vec3<Real>{px, py, pz};

            pxr_qs::generate_photon_update_momentum<Real, TableType, pxr::unit_system::SI>(
                chi, p_mom,
                get_rand<Real>::get(rand_gen),
                ref_table,
                e_phot);

            pdata.m_momentum(i, 0) = p_mom[0];
            pdata.m_momentum(i, 1) = p_mom[1];
            pdata.m_momentum(i, 2) = p_mom[2];
            photon_momentum(i, 0) = e_phot[0];
            photon_momentum(i, 1) = e_phot[1];
            photon_momentum(i, 2) = e_phot[2];

            rand_pool.free_state(rand_gen);
        });
    }

    return check_multi(photon_momentum, true, true) && check_multi(pdata.m_momentum, true, true);
}


template <typename Real>
void do_test(Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    auto particle_data = create_particles<Real>(
        how_many_particles,
        P_min, P_max, E_min, E_max, B_min, B_max, rand_pool);

    const auto dndt_table =
        generate_dndt_table<Real, KokkosVectorWrapper<Real>>(
            table_chi_min,
            table_chi_max,
            table_chi_size);

    const auto phot_em_table =
        generate_photon_emission_table<Real,KokkosVectorWrapper<Real>>(
            table_chi_min,
            table_chi_max,
            table_frac_min,
            table_chi_size,
            table_frac_size);

    const auto dndt_table_view = dndt_table.get_view();
    const auto phot_em_table_view = phot_em_table.get_view();

    const auto fill_opt_success =
        fill_opt_test<Real>(particle_data, how_many_repetitions, rand_pool);

    std::cout << ( fill_opt_success? "[ OK ]":"[ FAIL ]" )
        << std::endl;

    const auto evolve_opt_success =
        evolve_optical_depth<Real>(
            particle_data, dndt_table_view, dt_test, how_many_repetitions);

    std::cout << ( evolve_opt_success? "[ OK ]":"[ FAIL ]" )
        << std::endl;

    const auto phot_em_success =
        generate_photons<Real>(
            particle_data, phot_em_table_view, how_many_repetitions, rand_pool);

    std::cout << ( phot_em_success? "[ OK ]":"[ FAIL ]" )
        << std::endl;
}


int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    {
        Kokkos::Random_XorShift64_Pool<> rand_pool{random_seed};

        std::cout << "*** Kokkos example: begin ***" << std::endl;

        std::cout << "   --- Double precision test ---" << std::endl;
        do_test<double>(rand_pool);
        std::cout << "   --- END ---" << std::endl;

        std::cout << "   --- Single precision test ---" << std::endl;
        do_test<float>(rand_pool);
        std::cout << "   --- END ---" << std::endl;

        std::cout << "___ END ___" << std::endl;
    }
    Kokkos::finalize();
    return 0;
}
