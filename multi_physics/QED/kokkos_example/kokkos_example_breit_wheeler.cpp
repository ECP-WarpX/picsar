#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include <Kokkos_Random.hpp>

// PICSAR MULTIPHYSICS: BREIT WHEELER ENGINE
#define PXRMP_FORCE_INLINE
#define PXRMP_WITH_GPU
#define PXRMP_GPU_QUALIFIER KOKKOS_INLINE_FUNCTION

#include <picsar_qed/physics/chi_functions.hpp>
#include <picsar_qed/physics/gamma_functions.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_core.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp>
//__________________________________________________

#include <iostream>
#include <string>

//Some namespace aliases
namespace pxr =  picsar::multi_physics::phys;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_m =  picsar::multi_physics::math;
//__________________________________________________

//Some useful physical constants and functions
template<typename T = double>
constexpr T mec = static_cast<T>(pxr::electron_mass<> * pxr::light_speed<>);

template<typename T = double>
constexpr T mec2= static_cast<T>(pxr::electron_mass<> * pxr::light_speed<> * pxr::light_speed<>);

const double Es = pxr::schwinger_field<>;
//__________________________________________________

//Parameters of the test case
const unsigned int how_many_particles = 10'000'000;
const unsigned int how_many_repetitions = 1;
const double dt_test= 1e-18;
const double table_chi_min = 0.01;
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

template<typename Real>
struct ParticleData{
    static constexpr int num_components = 3;
    int num_particles = 0;
    Kokkos::View<Real * [num_components]> m_momentum;
    struct {
        Kokkos::View<Real *> Ex;
        Kokkos::View<Real *> Ey;
        Kokkos::View<Real *> Ez;
        Kokkos::View<Real *> Bx;
        Kokkos::View<Real *> By;
        Kokkos::View<Real *> Bz;
        Kokkos::View<Real *> opt;
    } m_fields;
};

template <typename Real>
class KokkosVectorWrapper : public Kokkos::vector<Real>
{
    public:
    template<typename... Args>
    KokkosVectorWrapper(Args&&... args) : Kokkos::vector<Real>(std::forward<Args>(args)...)
    {}

    void pxr_sync()
    {
        Kokkos::vector<Real>::on_device();
    }

    const Real* data() const
    {
        return Kokkos::vector<Real>::d_view.data();
    }

};

template<typename T>
std::string get_type_name()
{
    return typeid(T).name();
}

template <typename Real>
auto init_multi_comp_view_with_random_content(
    const std::string& label,
    const Real min_val, const Real max_val,
    const int N,
    Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    auto view = Kokkos::View<Real * [ParticleData<Real>::num_components]>{label, static_cast<unsigned int>(N)};

    Kokkos::fill_random(view, rand_pool, min_val, max_val);

    return view;
}

template <typename Real>
auto init_view_with_random_content(
    const std::string& label,
    const Real min_val, const Real max_val,
    const unsigned int N,
    Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    auto view = Kokkos::View<Real *>{label, N};
    Kokkos::fill_random(view, rand_pool, min_val, max_val);
    return view;
}

template<typename Real>
ParticleData<Real> create_particles(const int N, Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    ParticleData<Real> pdata;
    pdata.num_particles = N;

    pdata.m_momentum =
        init_multi_comp_view_with_random_content<Real>("mom", P_min, P_max, N, rand_pool);
    pdata.m_fields.Ex =
        init_view_with_random_content<Real>("Ex", E_min, E_max, N, rand_pool);
    pdata.m_fields.Ey =
        init_view_with_random_content<Real>("Ey", E_min, E_max, N, rand_pool);
    pdata.m_fields.Ez =
        init_view_with_random_content<Real>("Ez", E_min, E_max, N, rand_pool);
    pdata.m_fields.Bx =
        init_view_with_random_content<Real>("Bx", B_min, B_max, N, rand_pool);
    pdata.m_fields.By =
        init_view_with_random_content<Real>("By", B_min, B_max, N, rand_pool);
    pdata.m_fields.Bz =
        init_view_with_random_content<Real>("Bz", B_min, B_max, N, rand_pool);
    pdata.m_fields.opt =
        init_view_with_random_content<Real>("opt", 0.0, 0.0, N, rand_pool);
    return pdata;
}

template<typename Real>
void correct_low_momenta(ParticleData<Real>& pdata)
{
    const auto num_particles = pdata.num_particles;
    Kokkos::parallel_for("CorrectLowMomenta_"+get_type_name<Real>(),
            num_particles, KOKKOS_LAMBDA(int i){
            const auto px = pdata.m_momentum(i,0);
            const auto py = pdata.m_momentum(i,1);
            const auto pz = pdata.m_momentum(i,2);

            const auto gamma_gamma = pxr::compute_gamma_photon<Real>(px, py, pz);

            const auto bb = Real(2.1);

            if(gamma_gamma == Real(0.0) ){
                const auto cc = bb/std::sqrt(Real(3.0));
                pdata.m_momentum(i,0) = cc*mec<Real>;
                pdata.m_momentum(i,1) = cc*mec<Real>;
                pdata.m_momentum(i,2) = cc*mec<Real>;
            }
            else if (gamma_gamma < Real(2.0)){
                const auto cc = bb/gamma_gamma;
                pdata.m_momentum(i,0) *= cc;
                pdata.m_momentum(i,1) *= cc;
                pdata.m_momentum(i,2) *= cc;
            }
        });
}


template <typename Real, typename Vector>
auto generate_dndt_table(Real chi_min, Real chi_max, int chi_size)
{
    std::cout << "Preparing dndt table [" << typeid(Real).name() << ", " << chi_size <<"]...\n";
    std::cout.flush();

    pxr_bw::dndt_lookup_table_params<Real> bw_params{chi_min, chi_max, chi_size};

	auto table = pxr_bw::dndt_lookup_table<
        Real, Vector>{bw_params};

    table.generate();

    return table;
}

template <typename Real, typename Vector>
auto generate_pair_table(Real chi_min, Real chi_max, int chi_size, int frac_size)
{
    std::cout << "Preparing pair production table [" << typeid(Real).name() << ", " << chi_size << " x " << frac_size <<"]...\n";
    std::cout.flush();

    pxr_bw::pair_prod_lookup_table_params<Real> bw_params{
        chi_min, chi_max, chi_size, frac_size};

	auto table = pxr_bw::pair_prod_lookup_table<
        Real, Vector>{bw_params};

    table.template generate();

    return table;
}

template <typename Real>
bool check(Kokkos::View<Real *> field,
    const bool check_nan = true, const bool check_inf = false){

    int num = 0;

    if(check_nan && !check_inf){
        Kokkos::parallel_reduce("HowManyNaNs_"+get_type_name<Real>(),
        field.size(), KOKKOS_LAMBDA (const int& i, int& temp ) {
            temp += std::isnan(field(i));},
            num);
    }
    else if(!check_nan && check_inf){
        Kokkos::parallel_reduce("HowManyInfs_"+get_type_name<Real>(),
        field.size(), KOKKOS_LAMBDA (const int& i, int& temp ) {
            temp += std::isinf(field(i));},
            num);
    }
    else if(check_nan && check_inf){
        Kokkos::parallel_reduce("HowManyNaNsInfs_"+get_type_name<Real>(),
        field.size(), KOKKOS_LAMBDA (const int& i, int& temp ) {
            temp += std::isnan(field(i)) + std::isinf(field(i));},
            num);
    }

    return (num == 0);
}

template <typename Real>
bool check_multi(Kokkos::View<Real * [ParticleData<Real>::num_components]> vec,
    const bool check_nan = true, const bool check_inf = false)
{

    int num = 0;

    if(check_nan && !check_inf){
        Kokkos::parallel_reduce("HowManyNaNsMulti_"+get_type_name<Real>(),
        vec.extent(0), KOKKOS_LAMBDA (const int& i, int& temp ) {
            for (int cc = 0; cc < ParticleData<Real>::num_components; ++cc){
                temp += std::isnan(vec(i,cc));};},
            num);
    }
    else if(!check_nan && check_inf){
        Kokkos::parallel_reduce("HowManyInfsMulti_"+get_type_name<Real>(),
        vec.extent(0), KOKKOS_LAMBDA (const int& i, int& temp ) {
            for (int cc = 0; cc < ParticleData<Real>::num_components; ++cc){
                temp += std::isinf(vec(i,cc));};},
            num);
    }
    else if(check_nan && check_inf){
        Kokkos::parallel_reduce("HowManyNaNsInfsMulti_"+get_type_name<Real>(),
        vec.extent(0), KOKKOS_LAMBDA (const int& i, int& temp ) {
            for (int cc = 0; cc < ParticleData<Real>::num_components; ++cc){
                temp += std::isnan(vec(i,cc)) + std::isinf(vec(i,cc));};},
            num);
    }

    return (num == 0);
}

using GenType = Kokkos::Random_XorShift64_Pool<>::generator_type;

template <typename Real>
struct get_rand{
    KOKKOS_INLINE_FUNCTION
    static Real get(GenType& gen)
    {
        return gen.drand();
    }
};

template<>
struct get_rand<float>{

    KOKKOS_INLINE_FUNCTION
    static float get(GenType& gen)
    {
        return gen.frand();
    }
};

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
                pxr_bw::get_optical_depth<Real>(
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
            const auto ee = std::sqrt(px*px + py*py + pz*pz)*pxr::light_speed<Real>;
            const auto chi = pxr::chi_photon<Real, pxr::unit_system::SI>(
                px, py ,pz, ex, ey, ez, bx, by, bz);
            pxr_bw::evolve_optical_depth<Real, TableType>(
                ee, chi, dt, opt, ref_table);
        });
    }

    return check(pdata.m_fields.opt, true, false);
}


template <typename Real, typename TableType>
bool generate_pairs(
    ParticleData<Real>& pdata, const TableType& ref_table, const int repetitions,
    Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    const auto num_particles = pdata.num_particles;

    auto ele_momentum = init_multi_comp_view_with_random_content<Real>(
        "ele_momentum", 0.0, 0.0, num_particles, rand_pool);
    auto pos_momentum = init_multi_comp_view_with_random_content<Real>(
        "pos_momentum", 0.0, 0.0, num_particles, rand_pool);

    for(int rr = 0; rr < repetitions; ++rr){
        Kokkos::parallel_for("PairGen", num_particles, KOKKOS_LAMBDA(int i){
            const auto px = pdata.m_momentum(i,0);
            const auto py = pdata.m_momentum(i,1);
            const auto pz = pdata.m_momentum(i,2);
            const auto ex = pdata.m_fields.Ex(i);
            const auto ey = pdata.m_fields.Ey(i);
            const auto ez = pdata.m_fields.Ez(i);
            const auto bx = pdata.m_fields.Bx(i);
            const auto by = pdata.m_fields.By(i);
            const auto bz = pdata.m_fields.Bz(i);

            const auto chi = pxr::chi_photon<Real, pxr::unit_system::SI>(
                px, py ,pz, ex, ey, ez, bx, by, bz);

            auto rand_gen = rand_pool.get_state();

            auto e_mom = pxr_m::vec3<Real>{0,0,0};
            auto p_mom = pxr_m::vec3<Real>{0,0,0};

            pxr_bw::generate_breit_wheeler_pairs<Real, TableType, pxr::unit_system::SI>(
                chi, pxr_m::vec3<Real>{px, py, pz},
                get_rand<Real>::get(rand_gen),
                ref_table,
                e_mom, p_mom);

            ele_momentum(i, 0) = e_mom[0];
            ele_momentum(i, 1) = e_mom[1];
            ele_momentum(i, 2) = e_mom[2];
            pos_momentum(i, 0) = p_mom[0];
            pos_momentum(i, 1) = p_mom[1];
            pos_momentum(i, 2) = p_mom[2];

            rand_pool.free_state(rand_gen);
        });
    }

    return check_multi(ele_momentum, true, true) && check_multi(pos_momentum, true, true);
}


template <typename Real>
void do_test(Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    auto particle_data = create_particles<Real>(how_many_particles, rand_pool);
    correct_low_momenta(particle_data);

    const auto dndt_table =
        generate_dndt_table<Real, KokkosVectorWrapper<Real>>(
            table_chi_min,
            table_chi_max,
            table_chi_size);

    const auto pair_table =
        generate_pair_table<Real,KokkosVectorWrapper<Real>>(
            table_chi_min,
            table_chi_max,
            table_chi_size,
            table_frac_size);

    const auto dndt_table_view = dndt_table.get_view();
    const auto pair_table_view = pair_table.get_view();

    const auto fill_opt_success =
        fill_opt_test<Real>(particle_data, how_many_repetitions, rand_pool);

    std::cout << ( fill_opt_success? "[ OK ]":"[ FAIL ]" )
        << std::endl;

    const auto evolve_opt_success =
        evolve_optical_depth<Real>(
            particle_data, dndt_table_view, dt_test, how_many_repetitions);

    std::cout << ( evolve_opt_success? "[ OK ]":"[ FAIL ]" )
        << std::endl;

    const auto pair_prod_success =
        generate_pairs<Real>(
            particle_data, pair_table_view, how_many_repetitions, rand_pool);

    std::cout << ( pair_prod_success? "[ OK ]":"[ FAIL ]" )
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
