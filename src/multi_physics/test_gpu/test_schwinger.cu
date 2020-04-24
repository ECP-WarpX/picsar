#include <iostream>
#include <vector>
#include <array>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <random>

#include <cuda.h>
#include <curand.h>

#define PXRMP_GPU __host__ __device__
#define PXRMP_WITH_SI_UNITS
#include "../QED/src/physics/schwinger/schwinger_pair_engine_core.hpp"

//Alias for the picsar::multi_physics namespace
namespace pxr =  picsar::multi_physics::phys;
namespace pxrm =  picsar::multi_physics::math;

//How many test cases
const int test_size = 10'000'000;

//Tolerance for double precision calculations
const double double_tolerance = 5.0e-10;

//Tolerance for single precision calculations
const float float_tolerance = 5.0e-2;

//Templated tolerance
template <typename T>
T constexpr tolerance()
{
    if(std::is_same<T,float>::value)
        return float_tolerance;
    else
        return double_tolerance;
}

//*********************** SCHWINGER ENGINE: expected_pair_number ******************************
template <typename Real, pxr::unit_system UnitSystem>
__global__ void calc_expected_pair_number(
	int n,
    const Real* ex, const Real* ey, const Real* ez,
    const Real* bx, const Real* by, const Real* bz,
	Real volume, Real dt, Real* __restrict__ exp_pairs, Real ref = 1.0)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < n){
        exp_pairs[i] = pxr::schwinger::expected_pair_number<Real, UnitSystem>(
            ex[i], ey[i], ez[i], bx[i], by[i], bz[i], volume, dt, ref);
	}
}
//********************************************************************************************

template<typename Real, pxr::unit_system UnitSystem>
std::pair<bool, float> do_schwinger_test(
    const int size,
    const std::array<std::vector<double>,3 > & E_SI,
    const std::array<std::vector<double>,3 > & B_SI,
    const double t_volume,
    const double t_dt,
    const std::vector<double> & exp_sol,
    const double ref = 1.0)
{
    std::array<std::vector<Real>,3 > vE;
    std::array<std::vector<Real>,3 > vB;

    const auto factE = pxr::conv<
        pxr::quantity::E,
        pxr::unit_system::SI,
        UnitSystem,
        Real>::fact(1.0, ref);
    for (int i = 0; i < 3; ++i){
        vE[i].resize(size);
        std::transform(vE[i].begin(), vE[i].end(),
            E_SI[i].begin(), vE[i].begin(), [&](auto, auto a2){
                return static_cast<Real>(a2)*factE;
            });
    }

    const auto factB = pxr::conv<
        pxr::quantity::B,
        pxr::unit_system::SI,
        UnitSystem,
        Real>::fact(1.0, ref);
        for (int i = 0; i < 3; ++i){
            vB[i].resize(size);
            std::transform(vB[i].begin(), vB[i].end(),
                B_SI[i].begin(), vB[i].begin(), [&](auto, auto a2){
                    return static_cast<Real>(a2)*factB;
                });
        }

    const auto volume = t_volume*pxr::conv<pxr::quantity::volume, pxr::unit_system::SI,
        UnitSystem, Real>::fact(1.0, ref);
    const auto dt = t_dt*pxr::conv<pxr::quantity::time, pxr::unit_system::SI,
        UnitSystem, Real>::fact(1.0, ref);

    Real* d_ex;
    Real* d_ey;
    Real* d_ez;
    Real* d_bx;
    Real* d_by;
    Real* d_bz;
    Real* d_sol;
    const auto bytesize = size*sizeof(Real);
    cudaMalloc(&d_ex, bytesize);
    cudaMalloc(&d_ey, bytesize);
    cudaMalloc(&d_ez, bytesize);
    cudaMalloc(&d_bx, bytesize);
    cudaMalloc(&d_by, bytesize);
    cudaMalloc(&d_bz, bytesize);
    cudaMalloc(&d_sol, bytesize);
    cudaMemcpy(d_ex, vE[0].data(), bytesize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ey, vE[1].data(), bytesize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ez, vE[2].data(), bytesize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_bx, vB[0].data(), bytesize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_by, vB[1].data(), bytesize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_bz, vB[2].data(), bytesize, cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    calc_expected_pair_number<Real, UnitSystem>
        <<<(size+255)/256, 256>>>(
            size, d_ex, d_ey, d_ez, d_bx, d_by, d_bz,
            dt, volume, d_sol, ref);
    cudaEventRecord(stop);

    auto sol = std::vector<Real>(size);
    cudaMemcpy(sol.data(), d_sol, bytesize, cudaMemcpyDeviceToHost);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    bool is_ok = true;

    for(int i = 0; i < size; ++i)
    {
        if(fabs(static_cast<Real>(exp_sol[i])) <  tolerance<Real>()){
            if(fabs(static_cast<Real>(exp_sol[i]) - sol[i]) >  tolerance<Real>()){
                is_ok = false;
                break;
            }
        }
        else if(fabs(sol[i] - static_cast<Real>(exp_sol[i]))/static_cast<Real>(exp_sol[i]) > tolerance<Real>())
        {
            is_ok = false;
            break;
        }
    }

    cudaFree(d_ex);
    cudaFree(d_ey);
    cudaFree(d_ez);
    cudaFree(d_bx);
    cudaFree(d_by);
    cudaFree(d_bz);
    cudaFree(d_sol);

    return std::make_pair(is_ok, milliseconds);
}

int main()
{
    //Generate data SI

    const double dt = 1.0e-15;
    const double vol = 1.0e-27;

    auto E_SI = std::array<std::vector<double>,3 >
        {std::vector<double>(test_size),
        std::vector<double>(test_size),
        std::vector<double>(test_size)};
    auto B_SI = std::array<std::vector<double>,3 >
        {std::vector<double>(test_size),
        std::vector<double>(test_size),
        std::vector<double>(test_size)};
    auto sol = std::vector<double>(test_size);

    std::mt19937 rng{};
    std::uniform_real_distribution<double> unf{-10,10};

    for (auto& v : E_SI){
        std::generate(v.begin(), v.end(), [&](){
            return unf(rng)*pxr::schwinger_field<double>;
        });
    }

    for (auto& v : B_SI){
        std::generate(v.begin(), v.end(), [&](){
            return unf(rng)*pxr::schwinger_field<double>/pxr::light_speed<double>;
        });
    }


    for(int i = 0; i < test_size; ++i){
        sol[i] =
            pxr::schwinger::expected_pair_number<double, pxr::unit_system::SI>(
                E_SI[0][i], E_SI[1][i], E_SI[2][i],
                B_SI[0][i], B_SI[1][i], B_SI[2][i],
                vol, dt);
    }

    auto resprint = [](const std::string& tname,std::pair<bool, float> res){
        auto ok = res.first;
        auto time = res.second;
        std::cout << tname << ": ";
        if(ok){
            std::cout << time << " ms [ SUCCESS! ] \n";
        }
        else
            std::cout << time << " ms [ FAILED! ] \n";
    };

    resprint("SI, double", do_schwinger_test<double, pxr::unit_system::SI>(
        test_size, E_SI, B_SI, vol, dt, sol));
    resprint("SI, float", do_schwinger_test<float, pxr::unit_system::SI>(
        test_size, E_SI, B_SI, vol, dt, sol));

    const double ref_lambda = 800e-9;
    const double ref_omega = 2.0*pxrm::pi<double>*pxr::light_speed<double>/ref_lambda;

    resprint("OMEGA, double", do_schwinger_test<double, pxr::unit_system::norm_omega>(
                test_size, E_SI, B_SI, vol, dt, sol, ref_omega));
    resprint("OMEGA, float", do_schwinger_test<float, pxr::unit_system::norm_omega>(
                test_size, E_SI, B_SI, vol, dt, sol, ref_omega));

    resprint("LAMBDA, double", do_schwinger_test<double, pxr::unit_system::norm_lambda>(
                test_size, E_SI, B_SI, vol, dt, sol, ref_lambda));
    resprint("LAMBDA, float", do_schwinger_test<float, pxr::unit_system::norm_lambda>(
                test_size, E_SI, B_SI, vol, dt, sol, ref_lambda));

    resprint("HL, double", do_schwinger_test<double, pxr::unit_system::heaviside_lorentz>(
            test_size, E_SI, B_SI, vol, dt, sol));
    resprint("HL, float", do_schwinger_test<float, pxr::unit_system::heaviside_lorentz>(
            test_size, E_SI, B_SI, vol, dt, sol));


	return 0;
}
