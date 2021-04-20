
#include <Kokkos_Core.hpp>

#include <iostream>

// PICSAR MULTIPHYSICS: BREIT WHEELER ENGINE
#define PXRMP_WITH_GPU
#define PXRMP_GPU_QUALIFIER
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_core.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables.hpp>
#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp>
//__________________________________________________

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

const int HOW_MANY_PARTICLES = 10'000'000;
const int HOW_MANY_REPETITIONS = 10;

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

template <typename Real>
void do_test()
{
    particle_data<Real> particle_data;
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