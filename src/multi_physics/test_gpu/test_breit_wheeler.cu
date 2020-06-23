#include <iostream>
#include <vector>
#include <array>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <random>
#include <tuple>
#include <omp.h>

#include <cuda.h>
#include <curand.h>
#include <thrust/device_vector.h>

// PICSAR MULTIPHYSICS: BREIT WHEELER ENGINE
#define PXRMP_GPU __host__ __device__
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_core.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tables.hpp"
#include "../QED/src/physics/breit_wheeler/breit_wheeler_engine_tabulated_functions.hpp"
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
//const int test_size = 5;
const double dt_test = 1e-18;
const double dndt_chi_min = 0.01;
const double dndt_chi_max = 200.0;
const double dndt_frac_out_low = 0.1;
const double dndt_frac_out_high = 0.1;
const double max_intable_normalized_momentum = 1000.0;
const double max_intable_normalized_field = 0.02;
const int dndt_table_size_1 = 128;
const int dndt_table_size_2 = 256;
const double ref_lambda = 800e-9;
const double ref_omega = 2.0*pxr_m::pi<double>*pxr::light_speed<double>/ref_lambda;
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

    RealType* data()
    {
        return thrust::raw_pointer_cast(thrust::device_vector<RealType>::data());
    }
};
//__________________________________________________


//Particle data structure
template<typename Real>
struct pdata
{
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
template <typename Real, typename TableType, pxr::unit_system UnitSystem>
__global__ void bw_opt_depth_evol(
	int n,
	Real* p_data, Real dt, const TableType ref_table, Real ref = 1.0)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	const Real ref_mom = mec<Real>*pxr::conv<
        pxr::quantity::momentum,
        pxr::unit_system::SI,
        UnitSystem,
        Real>::fact(1.0, ref);
        
    const Real ref_en = mec2<Real>*pxr::conv<
        pxr::quantity::energy,
        pxr::unit_system::SI,
        UnitSystem,
        Real>::fact(1.0, ref);
	
	if (i < n){
	    const auto px = p_data[10*i+0];
	    const auto py = p_data[10*i+1];
	    const auto pz = p_data[10*i+2];
	    const auto ex = p_data[10*i+3];
	    const auto ey = p_data[10*i+4];
	    const auto ez = p_data[10*i+5];
	    const auto bx = p_data[10*i+6];
	    const auto by = p_data[10*i+7];
	    const auto bz = p_data[10*i+8];
	    auto opt = p_data[10*i+9];
	    	    
        const auto chi = pxr::chi_photon<Real, UnitSystem>(
            px, py, pz, ex, ey, ez, bx, by, bz, ref);
            
        const auto ppx = px/ref_mom;
	    const auto ppy = py/ref_mom;
	    const auto ppz = pz/ref_mom;        
        const auto gamma = sqrt(ppx*ppx + ppy*ppy + ppz*ppz);        
   
        pxr_bw::evolve_optical_depth<Real,TableType,UnitSystem>(
            ref_en*gamma, chi, dt, opt,
            ref_table, ref);
            
        //printf("GPU : %d %7.2e %7.2e %7.2e --> %7.2e \n", i, ref_en*gamma, chi, p_data[10*i+9], opt);
            
        p_data[10*i+9] = opt;
	}
}
//********************************************************************************************


template <typename RealType, pxr::unit_system UnitSystem>
std::vector<RealType>
transform_data(const std::vector<pdata<double>>& t_data, RealType ref_quantity = 1.0)
{
    const auto factP = pxr::conv<
        pxr::quantity::momentum,
        pxr::unit_system::SI,
        UnitSystem,
        RealType>::fact(1.0, ref_quantity);

    const auto factE = pxr::conv<
        pxr::quantity::E,
        pxr::unit_system::SI,
        UnitSystem,
        RealType>::fact(1.0, ref_quantity);
        
    const auto factB = pxr::conv<
        pxr::quantity::B,
        pxr::unit_system::SI,
        UnitSystem,
        RealType>::fact(1.0, ref_quantity);

    auto data_gpu = std::vector<RealType>(t_data.size()*10);

    for(int i = 0; i < t_data.size() ; ++i)
    {
        data_gpu[10*i+0] =  static_cast<RealType>(factP * t_data[i].px);
        data_gpu[10*i+1] =  static_cast<RealType>(factP * t_data[i].py);
        data_gpu[10*i+2] =  static_cast<RealType>(factP * t_data[i].pz);
        data_gpu[10*i+3] =  static_cast<RealType>(factE * t_data[i].ex);
        data_gpu[10*i+4] =  static_cast<RealType>(factE * t_data[i].ey);
        data_gpu[10*i+5] =  static_cast<RealType>(factE * t_data[i].ez);
        data_gpu[10*i+6] =  static_cast<RealType>(factB * t_data[i].bx);
        data_gpu[10*i+7] =  static_cast<RealType>(factB * t_data[i].by);
        data_gpu[10*i+8] =  static_cast<RealType>(factB * t_data[i].bz);
        data_gpu[10*i+9] = static_cast<RealType>(t_data[i].opt);
    }
    
    return data_gpu;
}


template <typename RealType, typename TableViewType, pxr::unit_system UnitSystem>
void do_dndt_test_units(
    TableViewType table,
    const std::vector<pdata<double>>& t_data, RealType dt, const std::vector<double>& sol,
    RealType ref_quantity = 1.0)
{
    auto how_many = t_data.size();
    auto data_gpu = transform_data<RealType, UnitSystem>(t_data);
    const auto bytesize = t_data.size()*10*sizeof(RealType);
    RealType* d_data;
    cudaMalloc(&d_data, bytesize);
    
    cudaMemcpy(d_data, data_gpu.data(), bytesize, cudaMemcpyHostToDevice);
    
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);    
    bw_opt_depth_evol<RealType, TableViewType, UnitSystem>
        <<<(how_many+255)/256, 256>>>(
            how_many, d_data, dt, table, ref_quantity);
    cudaEventRecord(stop);
            
    cudaMemcpy(data_gpu.data(), d_data, bytesize, cudaMemcpyDeviceToHost);
    
    int fail_count = 0;
    
    #pragma omp parallel for
    for(int i = 0; i < t_data.size(); ++i){
        const auto diff = data_gpu[10*i+9] - sol[i];
        const auto diff_norm = (sol[i] != 0.0 )?(diff):(diff/sol[i]); 
        if (diff_norm > tolerance<RealType>()){
            #pragma atomic
            fail_count++;
        }
    }

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << "      elapsed time : " << milliseconds << " ms ";
    if (fail_count == 0)
        std::cout << "[ OK! ]\n";
    else
        std::cout << "[ "<< fail_count*100.0/t_data.size() <<" % FAILS :-( ]\n"; 
    
    cudaFree(d_data);
}

template <typename RealType, typename TableViewType>
void do_dndt_test(
    TableViewType table, std::vector<pdata<double>> data, RealType dt,
    std::vector<double> sol, RealType rlambda, RealType romega)
{
    auto dt_l = dt * pxr::conv<
        pxr::quantity::time,
        pxr::unit_system::SI,
        pxr::unit_system::norm_lambda,
        RealType>::fact(1.0, rlambda);
        
    auto dt_o = dt * pxr::conv<
        pxr::quantity::time,
        pxr::unit_system::SI,
        pxr::unit_system::norm_omega,
        RealType>::fact(1.0, romega);
        
    auto dt_hl = dt * pxr::conv<
        pxr::quantity::time,
        pxr::unit_system::SI,
        pxr::unit_system::heaviside_lorentz,
        RealType>::fact();   
    
    std::cout << "      performing test in SI units:          "; std::cout.flush();
    do_dndt_test_units<RealType, TableViewType, pxr::unit_system::SI>(table, data, dt, sol);
    
    std::cout << "      performing test in norm_omega units:  "; std::cout.flush();
    do_dndt_test_units<RealType, TableViewType, pxr::unit_system::norm_omega>(table, data, dt_o, sol, romega);
    
    std::cout << "      performing test in norm_lambda units: "; std::cout.flush();
    do_dndt_test_units<RealType, TableViewType, pxr::unit_system::norm_lambda>(table, data, dt_l, sol, rlambda);
    
    std::cout << "      performing test in HL units:          "; std::cout.flush();
    do_dndt_test_units<RealType, TableViewType, pxr::unit_system::heaviside_lorentz>(table, data, dt_hl, sol);
}

template<typename TableType>
std::vector<double> dndt_calculate_solution(const std::vector<pdata<double>>& data, double dt, TableType table)
{
    std::vector<double> res(data.size());
       
    #pragma omp parallel for
    for(int i = 0; i < data.size(); ++i)
    {
        auto opt = data[i].opt;
        const auto chi = pxr::chi_photon<double, pxr::unit_system::SI>(
            data[i].px, data[i].py, data[i].pz, data[i].ex, data[i].ey, data[i].ez,  data[i].bx, data[i].by, data[i].bz);
            
        const double ppx = data[i].px/mec<>;
	    const double ppy = data[i].py/mec<>;
	    const double ppz = data[i].pz/mec<>;        
        const double en = mec2<>*sqrt(ppx*ppx + ppy*ppy + ppz*ppz);   
    
        pxr_bw::evolve_optical_depth<double, TableType, pxr::unit_system::SI>(
             en, chi, dt, opt, 
            table);
            
        //printf("CPU : %d %7.2e %7.2e %7.2e --> %7.2e \n", i, en, chi, data[i].opt, opt);
        res[i] = opt;
        
        
    } 
    return res;
}


auto generate_dndt_table_cpu(double chi_min, double chi_max, int table_size)
{
    std::cout << "Preparing dndt table for CPU..."; std::cout.flush();
    pxr_bw::dndt_lookup_table_params<double> bw_params{chi_min, chi_max, table_size};
	
	auto table = pxr_bw::dndt_lookup_table<
        double, std::vector<double>,
        pxr_bw::dndt_table_type::log,
		pxr_bw::dndt_table_out_policy::approx>{bw_params};
	
    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<double>(all_coords.size());
    
    #pragma omp parallel for
    for (int i = 0; i < all_vals.size(); ++i){
        all_vals[i] = pxr_bw::compute_T_function(all_coords[i]);
    }
    
    if(!table.set_all_vals(all_vals) )
        throw "Fail!";
        
    auto table_extrema = pxr_bw::dndt_lookup_table<
        double, std::vector<double>,
        pxr_bw::dndt_table_type::log,
		pxr_bw::dndt_table_out_policy::extrema>{bw_params};
		
    if(!table_extrema.set_all_vals(all_vals) )
        throw "Fail!";	
        
    std::cout << "done! \n"; std::cout.flush();
    return std::make_pair(table, table_extrema);
}


template <typename RealType>
auto generate_dndt_table_gpu(RealType chi_min, RealType chi_max, int table_size)
{
    std::cout << "Preparing dndt table [" << typeid(RealType).name() << ", " << table_size <<"]..."; std::cout.flush();
    pxr_bw::dndt_lookup_table_params<RealType> bw_params{chi_min, chi_max, table_size};
	
	auto table = pxr_bw::dndt_lookup_table<
        RealType, ThrustDeviceWrapper<RealType>,
        pxr_bw::dndt_table_type::log,
		pxr_bw::dndt_table_out_policy::approx>{bw_params};
	
    const auto all_coords = table.get_all_coordinates();
    auto all_vals = std::vector<RealType>(all_coords.size());
    
    #pragma omp parallel for
    for (int i = 0; i < all_vals.size(); ++i){
        all_vals[i] = pxr_bw::compute_T_function(all_coords[i]);
    }
    
    if(!table.set_all_vals(all_vals) )
        throw "Fail!";
        
    auto table_extrema = pxr_bw::dndt_lookup_table<
        RealType, ThrustDeviceWrapper<RealType>,
        pxr_bw::dndt_table_type::log,
		pxr_bw::dndt_table_out_policy::extrema>{bw_params};
		
    if(!table_extrema.set_all_vals(all_vals) )
        throw "Fail!";	
        
    std::cout << "done! \n"; std::cout.flush();
    return std::make_pair(table, table_extrema);
}


std::vector<pdata<double>> dndt_prepare_data(
    int how_many, double chi_min, double chi_max, double frac_out_low, double frac_out_high,
    double pscale, double fscale)
{
    std::cout << "Preparing data..."; std::cout.flush();    
    
    auto data = std::vector<pdata<double>>(how_many);
    
    std::mt19937 rng{};
    std::uniform_real_distribution<double> pp{-pscale*mec<>,pscale*mec<>};
    std::uniform_real_distribution<double> ee{ -Es*fscale, Es*fscale};
    std::uniform_real_distribution<double> bb{ -Bs*fscale, Bs*fscale};
    std::uniform_real_distribution<double> zero_one{0.0,1.0};
    auto exp_dist = std::exponential_distribution<double>(1.0);
    
    for(int i = 0; i < how_many; ++i){
        double chi = 0.0;
        double px,py,pz,ex,ey,ez,bx,by,bz,opt;
        do{
            px = pp(rng); py = pp(rng); pz = pp(rng);
            ex = ee(rng); ey = ee(rng); ez = ee(rng);
            bx = bb(rng); by = bb(rng); bz = bb(rng);
            chi = pxr::chi_photon<double, pxr::unit_system::SI>(
                px, py, pz, ex, ey, ez, bx, by, bz);
        } while(chi < chi_min || chi > chi_max);
        
        double cumulative_prob = zero_one(rng);
        
        if(cumulative_prob < frac_out_low){
            do{
                px *= 0.5; py *= 0.5; pz *= 0.5;
                chi = pxr::chi_photon<double, pxr::unit_system::SI>(
                    px, py, pz, ex, ey, ez, bx, by, bz);
            } while(chi >= chi_min);
        }
        else if(cumulative_prob < (frac_out_low+frac_out_high)){
            do{
                px = 2.0*px + mec<>; py = 2.0*py + mec<>; pz = 2.0*pz + mec<>;
                chi = pxr::chi_photon<double, pxr::unit_system::SI>(
                    px, py, pz, ex, ey, ez, bx, by, bz);
            } while(chi <= chi_max);        
        }
        opt = exp_dist(rng);
        data[i] = pdata<double>{px,py,pz,ex,ey,ez,bx,by,bz,opt};       
    }
    
    std::cout << "done!\n"; std::cout.flush();    
    
    return data;    
}

void do_dndt_test()
{
    std::cout << "--- dndt table ---\n";
     
    const auto data = dndt_prepare_data(
        test_size, dndt_chi_min, dndt_chi_max, dndt_frac_out_low, dndt_frac_out_high,
        max_intable_normalized_momentum, max_intable_normalized_field);

    auto tables_d1= generate_dndt_table_gpu<double>(dndt_chi_min, dndt_chi_max, dndt_table_size_1);
    auto tables_f1= generate_dndt_table_gpu<float>(dndt_chi_min, dndt_chi_max, dndt_table_size_1);
    auto tables_d2= generate_dndt_table_gpu<double>(dndt_chi_min, dndt_chi_max, dndt_table_size_2);
    auto tables_f2= generate_dndt_table_gpu<float>(dndt_chi_min, dndt_chi_max, dndt_table_size_2);
    
    auto dndt_d1_approx = tables_d1.first;
    auto dndt_d1_extrema = tables_d1.second;
    auto dndt_f1_approx = tables_f1.first;
    auto dndt_f1_extrema = tables_f1.second;
    auto dndt_d2_approx = tables_d2.first;
    auto dndt_d2_extrema = tables_d2.second;
    auto dndt_f2_approx = tables_f2.first;
    auto dndt_f2_extrema = tables_f2.second;
    
    auto v_dndt_d1_approx = dndt_d1_approx.get_view();
    auto v_dndt_d1_extrema= dndt_d1_extrema.get_view();
    auto v_dndt_f1_approx = dndt_f1_approx.get_view();
    auto v_dndt_f1_extrema= dndt_f1_extrema.get_view();
    auto v_dndt_d2_approx = dndt_d2_approx.get_view();
    auto v_dndt_d2_extrema= dndt_d2_extrema.get_view();
    auto v_dndt_f2_approx = dndt_f2_approx.get_view();
    auto v_dndt_f2_extrema= dndt_f2_extrema.get_view();
    
    auto tables_cpu_1 = generate_dndt_table_cpu(dndt_chi_min, dndt_chi_max, dndt_table_size_1);
    auto tables_cpu_2 = generate_dndt_table_cpu(dndt_chi_min, dndt_chi_max, dndt_table_size_2);
    auto dndt_dcpu_approx_1 = tables_cpu_1.first;
    auto dndt_dcpu_extrema_1 = tables_cpu_1.second;
    auto dndt_dcpu_approx_2 = tables_cpu_2.first;
    auto dndt_dcpu_extrema_2 = tables_cpu_2.second;
    std::cout << "Calculating solutions..."; std::cout.flush();
    const auto sol_approx_1 = dndt_calculate_solution(data, dt_test, dndt_dcpu_approx_1.get_view());
    const auto sol_extrema_1 = dndt_calculate_solution(data, dt_test, dndt_dcpu_extrema_1.get_view());
    const auto sol_approx_2 = dndt_calculate_solution(data, dt_test, dndt_dcpu_approx_2.get_view());
    const auto sol_extrema_2 = dndt_calculate_solution(data, dt_test, dndt_dcpu_extrema_2.get_view());
    std::cout << "done!\n"; std::cout.flush();    
    
    
    std::cout << "Performing test [double, table_size="<< dndt_table_size_1 <<", out_policy=approx]...\n"; std::cout.flush();
    do_dndt_test<double,decltype(v_dndt_d1_approx)>(v_dndt_d1_approx, data, dt_test, sol_approx_1, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();
    
    std::cout << "Performing test [float, table_size="<< dndt_table_size_1 <<", out_policy=approx]...\n"; std::cout.flush();
    do_dndt_test<float,decltype(v_dndt_f1_approx)>(v_dndt_f1_approx, data, dt_test, sol_approx_1, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();
    
    std::cout << "Performing test [double, table_size="<< dndt_table_size_1 <<", out_policy=extrema]...\n"; std::cout.flush();
    do_dndt_test<double,decltype(v_dndt_d1_extrema)>(v_dndt_d1_extrema, data, dt_test, sol_extrema_1, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();
    
    std::cout << "Performing test [float, table_size="<< dndt_table_size_1 <<", out_policy=extrema]...\n"; std::cout.flush();
    do_dndt_test<float,decltype(v_dndt_f1_extrema)>(v_dndt_f1_extrema, data, dt_test, sol_extrema_1, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();

    std::cout << "Performing test [double, table_size="<< dndt_table_size_2 <<", out_policy=approx]...\n"; std::cout.flush();
    do_dndt_test<double,decltype(v_dndt_d2_approx)>(v_dndt_d2_approx, data, dt_test, sol_approx_2, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();
 
    std::cout << "Performing test [float, table_size="<< dndt_table_size_2 <<", out_policy=approx]...\n"; std::cout.flush();
    do_dndt_test<float,decltype(v_dndt_f2_approx)>(v_dndt_f2_approx, data, dt_test, sol_approx_2, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();
    
    std::cout << "Performing test [double, table_size="<< dndt_table_size_2 <<", out_policy=extrema]...\n"; std::cout.flush();
    do_dndt_test<double,decltype(v_dndt_d2_extrema)>(v_dndt_d2_extrema, data, dt_test, sol_extrema_2, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();

    std::cout << "Performing test [float, table_size="<< dndt_table_size_2 <<", out_policy=extrema]...\n"; std::cout.flush();
    do_dndt_test<float,decltype(v_dndt_f2_extrema)>(v_dndt_f2_extrema, data, dt_test, sol_extrema_2, ref_lambda, ref_omega);
    std::cout << "...done! \n"; std::cout.flush();

    
    std::cout << "\n--- END ---\n";    
}

int main()
{
    std::cout << "*** Breit Wheeler engine GPU test *** \n \n";    
    do_dndt_test();    
    std::cout << "\n*** END ***\n";
	return 0;
}
