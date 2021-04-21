#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include <Kokkos_Random.hpp>

// PICSAR MULTIPHYSICS: BREIT WHEELER ENGINE
#define PXRMP_FORCE_INLINE
#define PXRMP_WITH_GPU
#define PXRMP_GPU_QUALIFIER KOKKOS_INLINE_FUNCTION
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

//Some useful physical constants
template<typename T = double>
constexpr T mec = static_cast<T>(pxr::electron_mass<> * pxr::light_speed<>);

template<typename T = double>
constexpr T mec2= static_cast<T>(pxr::electron_mass<> * pxr::light_speed<> * pxr::light_speed<>);

const double Es = pxr::schwinger_field<>;
const double Bs = pxr::schwinger_field<>/pxr::light_speed<>;
//__________________________________________________

//Parameters of the test case
const unsigned int how_many_particles = 10'000'000;
const unsigned int how_many_repetitions = 10;
const double dt_test = 1e-18;
const double table_chi_min = 0.01;
const double table_chi_max = 1000.0;
const int table_chi_size = 128;
const int table_frac_size = 128;
const double max_normalized_field = 0.02;
const double max_normalized_momentum = 1000.0;
const int random_seed = 22051988;
const double E_min = -0.01*Es;
const double E_max = 0.01*Es;
const double B_min = E_min/pxr::light_speed<>;
const double B_max = E_max/pxr::light_speed<>;
const double P_min = -mec<>;
const double P_max = mec<>;
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
    const int N,
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
bool check_not_nan(Kokkos::View<Real *> field){

   int num_nans = 0;
   Kokkos::parallel_reduce(
       "HowManyNaNs", field.size(),
        KOKKOS_LAMBDA (const int& i, int& temp ) {
           temp += std::isnan(field(i));},
        num_nans);


    return (num_nans == 0);
}

template <typename Real>
bool fill_opt_test(ParticleData<Real>& pdata, const int repetitions)
{
    for(int rr = 0; rr < repetitions; ++rr){
        const auto num_particles = pdata.num_particles;
        Kokkos::parallel_for("FillOpt", num_particles, KOKKOS_LAMBDA(int i){
            pdata.m_fields.opt(i) = pxr_bw::get_optical_depth<Real>(0.5);
        });
    }

    return check_not_nan(pdata.m_fields.opt);
}

template <typename Real>
void do_test(Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    auto particle_data = create_particles<Real>(how_many_particles, rand_pool);

    const auto dndt_table =
        generate_dndt_table<Real, Kokkos::vector<Real>>(
            table_chi_min,
            table_chi_max,
            table_chi_size);

    const auto pair_table =
        generate_pair_table<Real,Kokkos::vector<Real>>(
            table_chi_min,
            table_chi_max,
            table_chi_size,
            table_frac_size);

    const auto dndt_table_view = dndt_table.get_view();
    const auto pair_table_view = pair_table.get_view();

    bool fill_opt_success = fill_opt_test(particle_data, how_many_repetitions);

    std::cout << ( fill_opt_success? "[ OK ]":"[ FAIL ]" )
        << std::endl;
}


int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    {
        Kokkos::Random_XorShift64_Pool<> rand_pool{random_seed};

        std::cout << "*** Kokko example: begin ***" << std::endl;

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