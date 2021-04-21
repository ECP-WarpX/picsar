#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

// PICSAR MULTIPHYSICS: BREIT WHEELER ENGINE
#define PXRMP_FORCE_INLINE
#define PXRMP_WITH_GPU
#define PXRMP_GPU_QUALIFIER KOKKOS_INLINE_FUNCTION
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_core.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp>
//__________________________________________________

#include <iostream>

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
const int how_many_particles = 10'000'000;
const int how_many_repetitions = 10;
const double dt_test = 1e-18;
const double table_chi_min = 0.01;
const double table_chi_max = 1000.0;
const int table_chi_size = 256;
const int table_frac_size = 256;
const double max_normalized_field = 0.02;
const double max_normalized_momentum = 1000.0;
const int random_seed = 22051988;
//__________________________________________________


template<typename Real>
struct particle_data{
    static constexpr int num_components = 3;
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

template<typename Real>
particle_data<Real> create_particles(int how_many)
{
    particle_data<Real> pdata;

    pdata.m_momentum =
        Kokkos::View<Real * [particle_data<Real>::num_components]>{"mom", how_many};
    pdata.m_fields.Ex = Kokkos::View<Real *>{"Ex", how_many};
    pdata.m_fields.Ey = Kokkos::View<Real *>{"Ey", how_many};
    pdata.m_fields.Ez = Kokkos::View<Real *>{"Ez", how_many};
    pdata.m_fields.Bx = Kokkos::View<Real *>{"Bx", how_many};
    pdata.m_fields.By = Kokkos::View<Real *>{"By", how_many};
    pdata.m_fields.Bz = Kokkos::View<Real *>{"Bz", how_many};
    pdata.m_fields.opt = Kokkos::View<Real *>{"opt", how_many};
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
void do_test()
{
    particle_data<Real> particle_data;

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
}


int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);

    std::cout << "*** Kokko example: begin ***" << std::endl;

    std::cout << "   --- Double precision test ---" << std::endl;
    do_test<double>();
    std::cout << "   --- END ---" << std::endl;

    std::cout << "   --- Single precision test ---" << std::endl;
    do_test<float>();
    std::cout << "   --- END ---" << std::endl;

    std::cout << "___ END ___" << std::endl;

    Kokkos::finalize();
    return 0;
}