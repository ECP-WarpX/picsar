#ifndef __KOKKOS_EXAMPLE_COMMONS__
#define __KOKKOS_EXAMPLE_COMMONS__

//This file contains common functions and constants used by the two Kokkos examples

//Include Kokkos
#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include <Kokkos_Random.hpp>
//__________________________________________________

//Sets macros relevant for PICSAR QED according to the values of the corresponding Kokkos macros
#define PXRMP_FORCE_INLINE
#define PXRMP_WITH_GPU
#define PXRMP_GPU_QUALIFIER KOKKOS_INLINE_FUNCTION
//__________________________________________________

//Include PICSAR QED
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
#include <cstdlib>

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

/**
* Data structure to emulate particle data in a Particle-In-Cell code.
* Momenta are in a num_particlesx3 vector, while field components and optical depths
* are each in a separate vector of size num_particles.
*
* @tparam Real the floating point type to be used
*/
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

/**
* Thin Wrapper around Kokkos::vector<Real> to make it usable
* as the building block of the lookup tables.
*
* @tparam Real the floating point type to be used
*/
template <typename Real>
class KokkosVectorWrapper : public Kokkos::vector<Real>
{
    using KV = Kokkos::vector<Real>;
    public:

    /**
    * All the arguments passed to KokkosVectorWrapper constructor
    * are forwarded to Kokkos::vector constructor.
    *
    * @tparam Args the constructor arguments
    */
    template<typename... Args>
    KokkosVectorWrapper(Args&&... args) : KV(std::forward<Args>(args)...)
    {}

    /**
    * This function may be called in some places inside picsar_table classes.
    * It forces a deep copy of the CPU data to the GPU
    *
    */
    void pxr_sync()
    {
        Kokkos::deep_copy(KV::d_view, KV::h_view);
    }

    /**
    * This function returns a pointer to the raw data inside Kokkos::vector.
    * It is called to build table_views.
    *
    * @return a pointer to the raw vector data
    */
    const Real* data() const
    {
        return KV::d_view.data();
    }

};

/**
* An auxiliary function to call typeid(T).name()
*
* @tparam T a typename
* @return the name of T as a string
*/
template<typename T>
std::string get_type_name()
{
    return typeid(T).name();
}

/**
* Initializes a 2D Kokkos::View with N, ParticleData<Real>::num_components
* elements. The values inside the View are initialized randomly between min_val
* and max_val with a uniform distribution.
*
* @tparam Real the floating point type to be used
* @param[in] label the Kokkos::View label
* @param[in] min_val the minimum value to initialize the content of the Kokkos::View
* @param[in] max_val the maximum value to initialize the content of the Kokkos::View
* @param[in] N the number of particles
* @param[in,out] rand_pool a random pool
* @return the multi-component view
*/
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

/**
* Initializes a 1D Kokkos::View with N elements.
* The values inside the View are initialized randomly between min_val
* and max_val with a uniform distribution.
*
* @tparam Real the floating point type to be used
* @param[in] label the Kokkos::View label
* @param[in] min_val the minimum value to initialize the content of the Kokkos::View
* @param[in] max_val the maximum value to initialize the content of the Kokkos::View
* @param[in] N the number of particles
* @param[in,out] rand_pool a random pool
* @return the multi-component view
*/
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

/**
* Initializes "synthetic" particle data
* For each particle, the three components of momentum
* are randomly initialized with a uniform distribution between
* Pmin and Pmax, the three components of the electric field
* are initialized drawing from a uniform distribution between
* Emin and Emax, and the three components of the magnetic
* field are initialized drawing from a uniform distribution between
* Bmin and Bmax. The "optical depth" is initialized equal to zero
* for all the particles.
*
* @tparam Real the floating point type to be used
* @param N the number of particles
* @param[in] Pmin the minimum value of the momentum
* @param[in] Pmax the maximum value of the momentum
* @param[in] Emin the minimum value of the electric field
* @param[in] Emax the maximum value of the electric field
* @param[in] Bmin the minimum value of the magnetic field
* @param[in] Bmax the maximum value of the magnetic field
* @param[in,out] rand_pool a random pool
* @return the particle data
*/
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

/**
* Checks a 1D view for NaN and/or infs
*
* @tparam Real the floating point type to be used
* @param[in] field the view to be tested
* @param[in] check_nan whether to check for NaNs in field
* @param[in] check_inf whether to check for infinities in field
* @return true if checks pass, false otherwise
*/
template <typename Real>
bool check_nan_inf(Kokkos::View<Real *> field,
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

/**
* Checks a 2D view with N,ParticleData<Real>::num_components
* elements for NaN and/or infs
*
* @tparam Real the floating point type to be used
* @param[in] field the view to be tested
* @param[in] check_nan whether to check for NaNs in field
* @param[in] check_inf whether to check for infinities in field
* @return true if checks pass, false otherwise
*/
template <typename Real>
bool check_nan_inf_multi(Kokkos::View<Real * [ParticleData<Real>::num_components]> vec,
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

/**
* Calls either the drand() (double precision) of the frand() (single precision)
* method of the Kokkos random number generator, depending on the type Real.
* Prevents results exactly equal to 1.0
*
* @tparam Real the floating point type to be used
* @param[in,out] gen a random number generator
*/
template <typename Real>
KOKKOS_INLINE_FUNCTION
Real get_rand(GenType& gen)
{
    Real res = Real{1.0};

    if constexpr (std::is_same<Real,float>::value){
        while(res >= 1.0f)
            res = gen.frand();
    }
    else
    {
        while(res >= Real(1.0))
            res = gen.drand();
    }

    return res;
}

#endif //__KOKKOS_EXAMPLE_COMMONS__