#include <iostream>
#include <random>

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>

// PICSAR MULTIPHYSICS: BREIT WHEELER ENGINE
#define PXRMP_GPU __host__ __device__
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_core.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables_generator.hpp"
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

//Parameters of the test case
const int test_size = 20'000'000;
const double dt_test = 1e-18;
const double table_chi_min = 0.01;
const double table_chi_max = 500.0;
const int table_chi_size = 512;
const int table_frac_size = 512;
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
//__________________________________________________


//*********************** BREIT WHEELER ENGINE: optical depth evolution kernel ******************************
template <typename Real, typename TableType>
__global__ void bw_opt_depth_evol(
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
	    	    
        const auto chi = pxr::chi_photon<Real, pxr::unit_system::SI>(
            px, py, pz, ex, ey, ez, bx, by, bz);
            
        const auto ppx = px/mec<Real>;
	    const auto ppy = py/mec<Real>;
	    const auto ppz = pz/mec<Real>;        
        const auto en = (sqrt(1.0 + ppx*ppx + ppy*ppy + ppz*ppz))*mec2<Real>;        
   
        pxr_bw::evolve_optical_depth<Real,TableType>(
            en, chi, dt, opt,
            ref_table);            
        
        p_data[i].opt = opt;
	}
}
//********************************************************************************************

//*********************** BREIT WHEELER ENGINE: pair generation kernel ******************************
template <typename Real, typename TableType>
__global__ void bw_pair_gen(
	int n,
	const part<Real>* __restrict__ phot_data,
	const Real* __restrict__ random_numbers,
	part<Real>* __restrict__ ele_data,
	part<Real>* __restrict__ pos_data,
	const TableType ref_table)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < n){
	    const auto px = phot_data[i].px;
	    const auto py = phot_data[i].py;
	    const auto pz = phot_data[i].pz;
	    const auto ex = phot_data[i].ex;
	    const auto ey = phot_data[i].ey;
	    const auto ez = phot_data[i].ez;
	    const auto bx = phot_data[i].bx;
	    const auto by = phot_data[i].by;
	    const auto bz = phot_data[i].bz;
	    	    
        const auto chi = pxr::chi_photon<Real, pxr::unit_system::SI>(
            px, py, pz, ex, ey, ez, bx, by, bz);
                  
        const auto phot_mom = pxr_m::vec3<Real>{px, py, pz};
        auto ele_mom = pxr_m::vec3<Real>{};   
        auto pos_mom = pxr_m::vec3<Real>{};
        
        const auto unf_zero_one_minus_epsi = random_numbers[i];           
   
        pxr_bw::generate_breit_wheeler_pairs<Real,TableType, pxr::unit_system::SI>(
            chi, phot_mom, unf_zero_one_minus_epsi,
            ref_table,
            ele_mom, pos_mom);
        
        ele_data[i].x = phot_data[i].x;
        ele_data[i].y = phot_data[i].y;
        ele_data[i].z = phot_data[i].z;
        ele_data[i].px = ele_mom[0];
        ele_data[i].py = ele_mom[1];
        ele_data[i].pz = ele_mom[2];
        pos_data[i].ex = phot_data[i].x;
        pos_data[i].ey = phot_data[i].y;
        pos_data[i].ez = phot_data[i].z;
        pos_data[i].bx = pos_mom[0];
        pos_data[i].by = pos_mom[1];
        pos_data[i].bz = pos_mom[2]; 
	}
}
//********************************************************************************************

template <typename TableViewType>
void do_dndt_test(
    TableViewType table,
    const std::vector<part<double>>& t_data, double dt)
{
    auto data = t_data;
    auto how_many = data.size();
    const auto bytesize = t_data.size()*sizeof(part<double>);
    part<double>* d_data;
    cudaMalloc(&d_data, bytesize);
    
    cudaMemcpy(d_data, data.data(), bytesize, cudaMemcpyHostToDevice);
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);    
    bw_opt_depth_evol<double, TableViewType>
        <<<(how_many+255)/256, 256>>>(
            how_many, d_data, dt, table);
    cudaEventRecord(stop);
            
    cudaMemcpy(data.data(), d_data, bytesize, cudaMemcpyDeviceToHost);
    
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "      elapsed time : " << milliseconds << " ms \n";
    
    //hack
    if(data[0].px + data[how_many/2].py + data[how_many -1].pz != 0.0){
        std::cout.flush();
    }
      
    cudaFree(d_data);
}

template <typename TableViewType>
void do_pair_prod_test(
    TableViewType table,
    const std::vector<part<double>>& t_data, double dt)
{
    auto data = t_data;
    auto ele_data = t_data;
    auto pos_data = t_data;
    auto how_many = data.size();
    const auto bytesize = t_data.size()*sizeof(part<double>);
    part<double>* d_data;
    part<double>* d_ele_data;
    part<double>* d_pos_data;
    double* d_rand;
    cudaMalloc(&d_data, bytesize);
    cudaMalloc(&d_ele_data, bytesize);
    cudaMalloc(&d_pos_data, bytesize);
    cudaMalloc(&d_rand, sizeof(double)*how_many);
    
    cudaMemcpy(d_data, data.data(), bytesize, cudaMemcpyHostToDevice);
    
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, random_seed);
   
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    curandGenerateUniformDouble(gen, d_rand, how_many);  
    bw_pair_gen<double, TableViewType>
        <<<(how_many+255)/256, 256>>>(
            how_many, d_data, d_rand, d_ele_data, d_pos_data, table);
    cudaEventRecord(stop);
            
    cudaMemcpy(ele_data.data(), d_ele_data, bytesize, cudaMemcpyDeviceToHost);
    cudaMemcpy(pos_data.data(), d_pos_data, bytesize, cudaMemcpyDeviceToHost);
    
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "      elapsed time : " << milliseconds << " ms \n";
    
    //hack
    if(ele_data[0].px + pos_data[how_many/2].py + ele_data[how_many -1].pz != 0.0){
        std::cout.flush();
    }
      
    cudaFree(d_data);
    cudaFree(d_ele_data);
    cudaFree(d_pos_data);
    cudaFree(d_rand);
}

std::vector<part<double>>
prepare_data(
    int how_many, double pscale, double fscale)
{
    std::cout << "Preparing data..."; std::cout.flush();    
    
    auto data = std::vector<part<double>>(how_many);
    
    std::mt19937 rng{};
    std::uniform_real_distribution<double> pp{-pscale*mec<>,pscale*mec<>};
    std::uniform_real_distribution<double> ee{ -Es*fscale, Es*fscale};
    std::uniform_real_distribution<double> bb{ -Bs*fscale, Bs*fscale};
    std::uniform_real_distribution<double> pos{0.0,1.0};
    auto exp_dist = std::exponential_distribution<double>(1.0);
    
    for(int i = 0; i < how_many; ++i){
        data[i] = part<double>{
            pos(rng),pos(rng),pos(rng),
            pp(rng),pp(rng),pp(rng),
            ee(rng),ee(rng),ee(rng),
            bb(rng),bb(rng),bb(rng),
            exp_dist(rng)};               
    }
    
    std::cout << "done!\n"; std::cout.flush();    
    
    return data;    
}

template <typename RealType>
auto generate_dndt_table_gpu(RealType chi_min, RealType chi_max, int chi_size)
{
    std::cout << "Preparing dndt table [" << typeid(RealType).name() << ", " << chi_size <<"]...\n";
    std::cout.flush();
    
    pxr_bw::dndt_lookup_table_params<RealType> bw_params{chi_min, chi_max, chi_size};
	
	auto table = pxr_bw::dndt_lookup_table<
        RealType, ThrustDeviceWrapper<RealType>>{bw_params};
        
    table.template generate<pxr_bw::generation_policy::force_internal_double>();
	
    return table;
}

template <typename RealType>
auto generate_pair_table_gpu(RealType chi_min, RealType chi_max, int chi_size, int frac_size)
{
    std::cout << "Preparing pair production table [" << typeid(RealType).name() << ", " << chi_size << " x " << frac_size <<"]...\n";
    std::cout.flush();
    
    pxr_bw::pair_prod_lookup_table_params<RealType> bw_params{
        chi_min, chi_max, chi_size, frac_size};
	
	auto table = pxr_bw::pair_prod_lookup_table_logchi_linfrac<
        RealType, ThrustDeviceWrapper<RealType>>{bw_params};
        
    table.template generate<pxr_bw::generation_policy::force_internal_double>();
	
    return table;
}

void do_test()
{
    const auto data = prepare_data(test_size, max_normalized_momentum, max_normalized_field);
    const auto dndt_d_table = generate_dndt_table_gpu<double>(table_chi_min, table_chi_max, table_chi_size);
    const auto pair_d_table = generate_pair_table_gpu<double>(table_chi_min, table_chi_max, table_chi_size, table_frac_size);
    
    do_dndt_test(dndt_d_table.get_view(), data, dt_test);
    do_pair_prod_test(pair_d_table.get_view(), data, dt_test);
}

int main()
{
    std::cout << "*** Breit Wheeler engine GPU test *** \n \n";    
    do_test();    
    std::cout << "\n*** END ***\n";
	return 0;
}
