#ifndef __KOKKOS_EXAMPLE_COMMONS__
#define __KOKKOS_EXAMPLE_COMMONS__

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
#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_core.hpp>
#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp>
#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp>
//__________________________________________________

#include <iostream>
#include <string>

//Some namespace aliases
namespace pxr =  picsar::multi_physics::phys;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;
namespace pxr_m =  picsar::multi_physics::math;
//__________________________________________________


//Some useful physical constants and functions
template<typename T = double>
constexpr T mec = static_cast<T>(pxr::electron_mass<> * pxr::light_speed<>);

template<typename T = double>
constexpr T mec2= static_cast<T>(pxr::electron_mass<> * pxr::light_speed<> * pxr::light_speed<>);

const double Es = pxr::schwinger_field<>;
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
    using KV = Kokkos::vector<Real>;
    public:
    template<typename... Args>
    KokkosVectorWrapper(Args&&... args) : KV(std::forward<Args>(args)...)
    {}

    void pxr_sync()
    {
        Kokkos::deep_copy(KV::d_view, KV::h_view);
    }

    const Real* data() const
    {
        return KV::d_view.data();
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
ParticleData<Real> create_particles(const int N,
Real Pmin, Real Pmax, Real Emin, Real Emax, Real Bmin, Real Bmax,
Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
    ParticleData<Real> pdata;
    pdata.num_particles = N;

    pdata.m_momentum =
        init_multi_comp_view_with_random_content<Real>("mom", Pmin, Pmax, N, rand_pool);
    pdata.m_fields.Ex =
        init_view_with_random_content<Real>("Ex", Emin, Emax, N, rand_pool);
    pdata.m_fields.Ey =
        init_view_with_random_content<Real>("Ey", Emin, Emax, N, rand_pool);
    pdata.m_fields.Ez =
        init_view_with_random_content<Real>("Ez", Emin, Emax, N, rand_pool);
    pdata.m_fields.Bx =
        init_view_with_random_content<Real>("Bx", Bmin, Bmax, N, rand_pool);
    pdata.m_fields.By =
        init_view_with_random_content<Real>("By", Bmin, Bmax, N, rand_pool);
    pdata.m_fields.Bz =
        init_view_with_random_content<Real>("Bz", Bmin, Bmax, N, rand_pool);
    pdata.m_fields.opt =
        init_view_with_random_content<Real>("opt", 0.0, 0.0, N, rand_pool);
    return pdata;
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
        Real res = Real(1.0);
        while(res >= Real(1.0))
            res = gen.drand();
        return res;
    }
};

template<>
struct get_rand<float>{

    KOKKOS_INLINE_FUNCTION
    static float get(GenType& gen)
    {
        float res = 1.0f;
        while(res >= 1.0f)
            res = gen.frand();
        return res;
    }
};

#endif //__KOKKOS_EXAMPLE_COMMONS__