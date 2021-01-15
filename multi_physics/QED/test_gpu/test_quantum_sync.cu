/**
* This program tests the Quantum Synchrotron photon emission process on GPU.
* The test is performed in single and double precision. SI units
* are used for this test.
* Results are compared with those calculated on the CPU.
*/

#include <iostream>
#include <random>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>

// PICSAR MULTIPHYSICS: QUANTUM SYNCHROTRON ENGINE
#define PXRMP_GPU __host__ __device__
#include "../include/picsar_qed/physics/quantum_sync/quantum_sync_engine_core.hpp"
#include "../include/picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp"
#include "../include/picsar_qed/physics/quantum_sync/quantum_sync_engine_tables_generator.hpp"
//__________________________________________________

//Some namespace aliases
namespace pxr =  picsar::multi_physics::phys;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;
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
const int test_size = 100'000'000;
const double dt_test = 1e-18;
const double table_chi_min = 0.01;
const double table_chi_max = 1000.0;
const double table_frac_min = 1.0e-12;
const int table_chi_size = 256;
const int table_frac_size = 256;
const double max_normalized_field = 0.02;
const double max_normalized_momentum = 1000.0;
const int random_seed = 22051988;
//__________________________________________________

//Templated tolerance for double and single precision
const double double_tolerance = 5.0e-6;
const float float_tolerance = 5.0e-3;
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//Templated small number for double and single precision
const double double_small = 1.0e-10;
const float float_small = 1.0e-4;
template <typename T>
T constexpr small()
{
    if(std::is_same<T,float>::value)
        return float_small;
    else
        return double_small;
}
//__________________________________________________

//A thin wrapper around thrust::device_vector
template<typename RealType>
class ThrustDeviceWrapper : public thrust::device_vector<RealType>
{
    public:
    template<typename... Args>
    ThrustDeviceWrapper(Args&&... args) : thrust::device_vector<RealType>(std::forward<Args>(args)...)
    {}

    const RealType* data() const
    {
        return thrust::raw_pointer_cast(thrust::device_vector<RealType>::data());
    }
};
//__________________________________________________

//Particle data structure
template<typename Real>
struct part
{
    Real x;
    Real y;
    Real z;
    Real px;
    Real py;
    Real pz;
    Real ex;
    Real ey;
    Real ez;
    Real bx;
    Real by;
    Real bz;
    Real opt;
};

template<typename TIN, typename TOUT>
part<TOUT>
part_transform (const part<TIN>& in)
{
    part<TOUT> out;
    out.x = in.x ; out.y = in.y ; out.z = in.z ;
    out.px = in.px ; out.py = in.py ; out.pz = in.pz ;
    out.ex = in.ex ; out.ey = in.ey ; out.ez = in.ez ;
    out.bx = in.bx ; out.by = in.by ; out.bz = in.bz ;
    out.opt = in.opt ;
    return out;
}
//__________________________________________________

//*********************** QUANTUM SYNCHROTRON ENGINE: optical depth evolution kernel ******************************
template <typename Real, typename TableType>
__global__ void qs_opt_depth_evol(
	int n,
	part<Real>* p_data, Real dt, const TableType ref_table)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < n){
            const auto px = p_data[i].px;
            const auto py = p_data[i].py;
            const auto pz = p_data[i].pz;
            const auto ex = p_data[i].ex;
            const auto ey = p_data[i].ey;
            const auto ez = p_data[i].ez;
            const auto bx = p_data[i].bx;
            const auto by = p_data[i].by;
            const auto bz = p_data[i].bz;
            auto opt = p_data[i].opt;

            const auto chi = pxr::chi_ele_pos<Real, pxr::unit_system::SI>(
                px, py, pz, ex, ey, ez, bx, by, bz);

            const auto ppx = px/mec<Real>;
            const auto ppy = py/mec<Real>;
            const auto ppz = pz/mec<Real>;
            const auto en = (sqrt(static_cast<Real>(1.0) + ppx*ppx + ppy*ppy + ppz*ppz))*mec2<Real>;

            pxr_qs::evolve_optical_depth<Real,TableType>(
                en, chi, dt, opt,
                ref_table);

            p_data[i].opt = opt;
	}
}

template <typename Real, typename TableType>
void cpu_check_depth_evol_kernel(part<Real>& p_data, Real dt, const TableType ref_table)
{
    const auto px = p_data.px;
	const auto py = p_data.py;
	const auto pz = p_data.pz;
	const auto ex = p_data.ex;
	const auto ey = p_data.ey;
	const auto ez = p_data.ez;
	const auto bx = p_data.bx;
	const auto by = p_data.by;
	const auto bz = p_data.bz;
	auto opt = p_data.opt;

    const auto chi = pxr::chi_ele_pos<Real, pxr::unit_system::SI>(
            px, py, pz, ex, ey, ez, bx, by, bz);

    const auto ppx = px/mec<Real>;
	const auto ppy = py/mec<Real>;
	const auto ppz = pz/mec<Real>;
    const auto en = (sqrt(static_cast<Real>(1.0) + ppx*ppx + ppy*ppy + ppz*ppz))*mec2<Real>;

     pxr_qs::evolve_optical_depth<Real,TableType>(
        en, chi, dt, opt,
        ref_table);

    p_data.opt = opt;
}

//********************************************************************************************

//*********************** QUANTUM SYNCHROTRON ENGINE: photon emission kernel ******************************
template <typename Real, typename TableType>
__global__ void qs_phot_em(
	int n,
	part<Real>* __restrict__ p_data,
	const Real* __restrict__ random_numbers,
	part<Real>* __restrict__ phot_data,
	const TableType ref_table)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < n){
            const auto px = p_data[i].px;
            const auto py = p_data[i].py;
            const auto pz = p_data[i].pz;
            const auto ex = p_data[i].ex;
            const auto ey = p_data[i].ey;
            const auto ez = p_data[i].ez;
            const auto bx = p_data[i].bx;
            const auto by = p_data[i].by;
            const auto bz = p_data[i].bz;

            const auto chi = pxr::chi_ele_pos<Real, pxr::unit_system::SI>(
                px, py, pz, ex, ey, ez, bx, by, bz);

            auto p_mom = pxr_m::vec3<Real>{px, py, pz};
            auto phot_mom = pxr_m::vec3<Real>{};

            const auto unf_zero_one_minus_epsi = random_numbers[i];

            pxr_qs::generate_photon_update_momentum<Real,TableType, pxr::unit_system::SI>(
                chi, p_mom, unf_zero_one_minus_epsi,
                ref_table,
                phot_mom);

            p_data[i].px = p_mom[0];
            p_data[i].py = p_mom[1];
            p_data[i].pz = p_mom[2];
            phot_data[i].x = p_data[i].x;
            phot_data[i].y = p_data[i].y;
            phot_data[i].z = p_data[i].z;
            phot_data[i].px = phot_mom[0];
            phot_data[i].py = phot_mom[1];
            phot_data[i].pz = phot_mom[2];
	}
}

template <typename Real, typename TableType>
void cpu_check_phot_em_kernel(
    part<Real>& p_data,
    Real random_number,
    part<Real>& phot_data,
    const TableType ref_table)
{
    const auto px = p_data.px;
	const auto py = p_data.py;
	const auto pz = p_data.pz;
	const auto ex = p_data.ex;
	const auto ey = p_data.ey;
	const auto ez = p_data.ez;
	const auto bx = p_data.bx;
	const auto by = p_data.by;
	const auto bz = p_data.bz;

    const auto chi = pxr::chi_ele_pos<Real, pxr::unit_system::SI>(
        px, py, pz, ex, ey, ez, bx, by, bz);

    auto p_mom = pxr_m::vec3<Real>{px, py, pz};
    auto phot_mom = pxr_m::vec3<Real>{};

    const auto unf_zero_one_minus_epsi = random_number;

    pxr_qs::generate_photon_update_momentum<Real,TableType, pxr::unit_system::SI>(
        chi, p_mom, unf_zero_one_minus_epsi,
        ref_table,
        phot_mom);

    p_data.px = p_mom[0];
    p_data.py = p_mom[1];
    p_data.pz = p_mom[2];
    phot_data.x = p_data.x;
    phot_data.y = p_data.y;
    phot_data.z = p_data.z;
    phot_data.px = phot_mom[0];
    phot_data.py = phot_mom[1];
    phot_data.pz = phot_mom[2];
}

//********************************************************************************************

template <typename RealType, typename TableViewType>
void do_dndt_test(
    TableViewType table,
    const std::vector<part<double>>& t_data,
    double dt, TableViewType cpu_table)
{
    auto data = std::vector<part<RealType>>(t_data.size());
    std::transform(t_data.begin(), t_data.end(), data.begin(),
            part_transform<double, RealType>);

    auto how_many = data.size();
    const auto bytesize = t_data.size()*sizeof(part<RealType>);
    part<RealType>* d_data;
    cudaMalloc(&d_data, bytesize);

    cudaMemcpy(d_data, data.data(), bytesize, cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    qs_opt_depth_evol<RealType, TableViewType>
        <<<(how_many+511)/512, 512>>>(
            how_many, d_data, dt, table);
    cudaEventRecord(stop);

    cudaMemcpy(data.data(), d_data, bytesize, cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    float milliseconds = 0;

    bool is_ok = true;
    #pragma omp parallel for
    for(int i = 0; i < t_data.size(); ++i)
    {
        if(!is_ok) continue;

        const auto res_opt = data[i].opt;
        auto original_data = part_transform<double, RealType>(t_data[i]);
        cpu_check_depth_evol_kernel(original_data, static_cast<RealType>(dt), cpu_table);
        const auto exp_opt = original_data.opt;

        if(exp_opt < small<RealType>()){
            if(fabs(exp_opt - res_opt) >  small<RealType>()){
                is_ok = false;
                #pragma omp critical
                {
                    std::cout << exp_opt << " " << res_opt << std::endl;
                }
            }
        }
        else if(fabs((exp_opt - res_opt)/exp_opt) > tolerance<RealType>())
        {
            is_ok = false;
            #pragma omp critical
            {
                std::cout << exp_opt << " " << res_opt << std::endl;
            }
        }
    }

    cudaEventElapsedTime(&milliseconds, start, stop);
    if (is_ok)
        std::cout << "  [OK]  ";
    else
        std::cout << " [FAIL] ";
    std::cout << "elapsed time : " << milliseconds << " ms \n";

    cudaFree(d_data);
}

template <typename RealType, typename TableViewType>
void do_phot_em_test(
    TableViewType table,
    const std::vector<part<double>>& t_data, double dt, TableViewType cpu_table)
{
    auto data = std::vector<part<RealType>>(t_data.size());
    std::transform(t_data.begin(), t_data.end(), data.begin(),
            part_transform<double, RealType>);

    auto phot_data = t_data;
    auto how_many = data.size();
    auto rand_nums = std::vector<RealType>(data.size());
    const auto bytesize = t_data.size()*sizeof(part<RealType>);
    part<RealType>* d_data;
    part<RealType>* d_phot_data;
    RealType* d_rand;
    cudaMalloc(&d_data, bytesize);
    cudaMalloc(&d_phot_data, bytesize);
    cudaMalloc(&d_rand, sizeof(RealType)*how_many);

    cudaMemcpy(d_data, data.data(), bytesize, cudaMemcpyHostToDevice);

    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, random_seed);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    if(std::is_same<RealType, double>::value)
        curandGenerateUniformDouble(gen, reinterpret_cast<double*>(d_rand), how_many);
    else
        curandGenerateUniform(gen, reinterpret_cast<float*>(d_rand), how_many);

    qs_phot_em<RealType, TableViewType>
        <<<(how_many+511)/512, 512>>>(
            how_many, d_data, d_rand, d_phot_data, table);
    cudaEventRecord(stop);

    cudaMemcpy(phot_data.data(), d_phot_data, bytesize, cudaMemcpyDeviceToHost);
    cudaMemcpy(rand_nums.data(), d_rand, sizeof(RealType)*rand_nums.size(), cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    float milliseconds = 0;

    bool is_ok = true;
    #pragma omp parallel for
    for(int i = 0; i < t_data.size(); ++i)
    {
        if(!is_ok) continue;

        auto res_mom_phot = pxr_m::vec3<RealType>{
            phot_data[i].px,
            phot_data[i].py,
            phot_data[i].pz};
        auto original_data = part_transform<double, RealType>(t_data[i]);
        part<RealType> exp_phot_data;
        cpu_check_phot_em_kernel(original_data, rand_nums[i], exp_phot_data, cpu_table);
        auto exp_mom_phot = pxr_m::vec3<RealType>{
            exp_phot_data.px,
            exp_phot_data.py,
            exp_phot_data.pz};

        using namespace pxr_m;

        auto diff = res_mom_phot-exp_mom_phot;
        const auto dd = pxr_m::norm<RealType>(diff);

        if(fabs(dd) >  tolerance<RealType>()){
            is_ok = false;
        }
    }

    cudaEventElapsedTime(&milliseconds, start, stop);
    if (is_ok)
        std::cout << "  [OK]  ";
    else
        std::cout << " [FAIL] ";

    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "      elapsed time : " << milliseconds << " ms \n";

    cudaFree(d_data);
    cudaFree(d_phot_data);
    cudaFree(d_rand);
}


template<typename RealType>
std::vector<part<RealType>>
prepare_data(
    int how_many, RealType pscale, RealType fscale)
{
    std::cout << "Preparing data..."; std::cout.flush();

    auto data = std::vector<part<RealType>>(how_many);

    std::mt19937 rng{};
    std::uniform_real_distribution<RealType> pp{RealType(-pscale*mec<>),RealType(pscale*mec<>)};
    std::uniform_real_distribution<RealType> ee{RealType(-Es*fscale), RealType(Es*fscale)};
    std::uniform_real_distribution<RealType> bb{RealType(-Bs*fscale), RealType(Bs*fscale)};
    std::uniform_real_distribution<RealType> pos{0.0,1.0};
    auto exp_dist = std::exponential_distribution<RealType>(1.0);

    for(int i = 0; i < how_many; ++i){
        data[i] = part<RealType>{
            pos(rng),pos(rng),pos(rng),
            pp(rng),pp(rng),pp(rng),
            ee(rng),ee(rng),ee(rng),
            bb(rng),bb(rng),bb(rng),
            exp_dist(rng)};
    }

    std::cout << "done!\n"; std::cout.flush();

    return data;
}

template <typename RealType, typename VectorType>
auto generate_dndt_table(RealType chi_min, RealType chi_max, int chi_size)
{
    std::cout << "Preparing dndt table [" << typeid(RealType).name() << ", " << chi_size <<"]...\n";
    std::cout.flush();

    pxr_qs::dndt_lookup_table_params<RealType> qs_params{chi_min, chi_max, chi_size};

	auto table = pxr_qs::dndt_lookup_table<
        RealType, VectorType>{qs_params};

    table.generate();

    return table;
}

template <typename RealType, typename VectorType>
auto generate_phot_em_table(RealType chi_min, RealType chi_max, RealType frac_min, int chi_size, int frac_size)
{
    std::cout << "Preparing photon emission table [" << typeid(RealType).name() << ", " << chi_size << " x " << frac_size <<"]...\n";
    std::cout.flush();

    pxr_qs::photon_emission_lookup_table_params<RealType> qs_params{
        chi_min, chi_max, frac_min, chi_size, frac_size};

	auto table = pxr_qs::photon_emission_lookup_table<
        RealType, VectorType>{qs_params};

    table.template generate();

    return table;
}

template<typename RealType>
void do_test(const std::vector<part<double>>& data)
{
    std::cout << "Generating tables...\n"; std::cout.flush();
    const auto dndt_table =
        generate_dndt_table<RealType, ThrustDeviceWrapper<RealType>>(
            table_chi_min,
            table_chi_max,
            table_chi_size);

    const auto dndt_table_cpu =
        generate_dndt_table<RealType,std::vector<RealType>>(
            table_chi_min,
            table_chi_max,
            table_chi_size);

    const auto photem_table =
        generate_phot_em_table<RealType,ThrustDeviceWrapper<RealType>>(
            table_chi_min,
            table_chi_max,
            table_frac_min,
            table_chi_size,
            table_frac_size);

    const auto photem_table_cpu =
        generate_phot_em_table<RealType,std::vector<RealType>>(
            table_chi_min,
            table_chi_max,
            table_frac_min,
            table_chi_size,
            table_frac_size);

    std::cout << "done!\n"; std::cout.flush();

    std::cout << "Performing tests...\n"; std::cout.flush();
    do_dndt_test<RealType,typeof(dndt_table.get_view())>(dndt_table.get_view(), data, dt_test, dndt_table_cpu.get_view());
    do_phot_em_test<RealType,typeof(photem_table.get_view())>(photem_table.get_view(), data, dt_test, photem_table_cpu.get_view());
    std::cout << "done!\n"; std::cout.flush();
}

int main()
{
    std::cout << "*** Quantum Synchrotron engine GPU test *** \n \n"; std::cout.flush();
    const auto data = prepare_data<double>(
        test_size,
        max_normalized_momentum,
        max_normalized_field);
    std::cout << "*** Performing test in double precision ***\n"; std::cout.flush();
    do_test<double>(data);
    std::cout << "*** Performing test in single precision ***\n"; std::cout.flush();
    do_test<float>(data);
    std::cout << "\n*** END ***\n";
	return 0;
}
