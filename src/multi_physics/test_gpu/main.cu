#include <iostream>
#include <vector>

#include <cstdio>

#include <cuda.h>
#include <curand.h>

#define PXRMP_GPU __host__ __device__
#define PXRMP_WITH_SI_UNITS
#include "../QED/src/breit_wheeler_engine.hpp"
#include "../QED/src/rng_wrapper.hpp"

//Alias for the picsar::multi_physics namespace
namespace pxrmp =  picsar::multi_physics;

//Seed of the random number generator
const size_t seed = 83871734;

//A lot of particles!
const size_t N = 4000000;

//Sampling parameter for BW pair generation
const size_t sampling = 4;

//How many times should we repeat?
const int repeat = 1;

//Physical constants
const double las_wavlngth = 800.0 * pxrmp::si_nanometer;
const double me_c = pxrmp::electron_mass * pxrmp::light_speed;
const double eref = 2.0*pxrmp::pi*pxrmp::electron_mass*pxrmp::light_speed*pxrmp::light_speed/
            (las_wavlngth*pxrmp::elementary_charge);
const double bref = eref/pxrmp::light_speed;
//____________________________

//Const for momentum initialization
const double mom_coeff = 1000.0*me_c;

//Const for field initialization
const double efl_coeff = 1000.0*eref;
const double bfl_coeff = 1000.0*bref;

//An empty class to be used to call static BW functions
class dummy{};

//Lambda is not needed if SI units are used
const double default_lambda = 1.0;

//timestep to be used in this "simulation"
const double timestep = 1.0e-16;

//Helper function
bool does_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}



//***********************	BW ENGINE: get_optical_depth	******************************
//GPU kernel to initialize an array of optical depths given an
//array of random numbers [0,1)
__global__
void init_opt_depth(int n, double* opt, double* rand_nums)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < n){
		opt[i] = pxrmp::breit_wheeler_engine<double, dummy>::
			internal_get_optical_depth(1.0-rand_nums[i]);
	}
}
//********************************************************************************************



//GPU kernel to initialize fields & momenta randomly
__global__
void init_mom_fields(int n, double* px, double* py, double* pz,
	double*ex, double* ey, double* ez, double* bx, double* by, double* bz, double* ww, double* rand_nums)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < n){
		px[i] = mom_coeff*(0.5 - rand_nums[0*n + i]);
		py[i] = mom_coeff*(0.5 - rand_nums[1*n + i]);
		pz[i] = mom_coeff*(0.5 - rand_nums[2*n + i]);
		ex[i] = efl_coeff*(0.5 - rand_nums[3*n + i]);
		ey[i] = efl_coeff*(0.5 - rand_nums[4*n + i]);
		ez[i] = efl_coeff*(0.5 - rand_nums[5*n + i]);
		bx[i] = bfl_coeff*(0.5 - rand_nums[6*n + i]);
		by[i] = bfl_coeff*(0.5 - rand_nums[7*n + i]);
		bz[i] = bfl_coeff*(0.5 - rand_nums[8*n + i]);
		ww[i] = 1.0; // All weights are set to 1
	}
}




//*********************** BW ENGINE: evolve_opt_depth_and_determine_event ******************************
//GPU kernel to test internal_evolve_opt_depth_and_determine_event
__global__
void test_internal_evolve_opt_depth_and_determine_event(
	int n, double* px, double* py, double* pz, double* ex, double* ey, double* ez, double* bx, double* by, double* bz,
	double dt, double* opt, bool* has_event_happened, double* event_dt,
	size_t tab_how_many, double* coords, double* data, pxrmp::breit_wheeler_engine_ctrl<double>* bw_ctrl
)
{
	//Regenerate the lookuptable on GPU
	//This constructor does NOT allocate new memory: it manages existing pointers.
	pxrmp::lookup_1d<double> TTfunctab{tab_how_many, coords, data};

	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < n){
		pxrmp::breit_wheeler_engine<double, dummy>::
		internal_evolve_opt_depth_and_determine_event(
         		px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i],
         		dt, opt[i], has_event_happened[i], event_dt[i],
         		default_lambda, TTfunctab, *bw_ctrl);

	}
}
//********************************************************************************************


//GPU kernel to set all the has_event_happened to TRUE
__global__
void set_all_events_true(int n, bool* has_event_happened)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < n)
		has_event_happened[i] = true;
}

//*********************** BW ENGINE: generate_breit_wheeler_pairs ******************************
//GPU kernel to test internal_generate_breit_wheeler_pairs
__global__
void test_internal_generate_breit_wheeler_pairs(int n, bool* has_event_happened,
	double* px, double* py, double* pz, double* ex, double* ey, double* ez, double* bx, double* by, double* bz,
	double* weight, size_t sampling,
	double* e_px, double* e_py, double* e_pz, double* p_px, double* p_py, double* p_pz, double* e_weight, double* p_weight,
	size_t tab_how_many_1, double* coords_1, size_t tab_how_many_2, double* coords_2, double* data,
	pxrmp::breit_wheeler_engine_ctrl<double>* bw_ctrl, double* rand_num)
{
	pxrmp::lookup_2d<double> cum_prob_tab{tab_how_many_1, coords_1, tab_how_many_2, coords_2, data};

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int si = sampling*i;

	if (i < n && has_event_happened[i]){
		pxrmp::breit_wheeler_engine<double, dummy>::
		internal_generate_breit_wheeler_pairs(
			px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i], weight[i], sampling,
			&e_px[si], &e_py[si], &e_pz[si],
			&p_px[si], &p_py[si], &p_pz[si],
			&e_weight[si], &p_weight[si],
			default_lambda, cum_prob_tab, *bw_ctrl, &rand_num[si]);
	}
}

//********************************************************************************************

int main()
{
	//Seed will be used only with cuRand
	size_t useless_seed = 22051988;

	//Lambda is not used with SI units
	double useless_lambda = 1.0;

	//Change default table parameters in order to speed up the calculations
	pxrmp::breit_wheeler_engine_ctrl<double> bw_ctrl;
	bw_ctrl.chi_phot_tdndt_how_many = 200;
	bw_ctrl.chi_phot_tpair_how_many = 3;
	bw_ctrl.chi_frac_tpair_how_many = 3;

	//Initialize the BW engine
	auto bw_engine =
		pxrmp::breit_wheeler_engine<double, pxrmp::stl_rng_wrapper>
		{std::move(pxrmp::stl_rng_wrapper{useless_seed}), useless_lambda, bw_ctrl};

	//Initialize the lookup tables
   	//Generates tables if they do not exist
   	if(!does_file_exist("tdndt.bin")){
       		bw_engine.compute_dN_dt_lookup_table(&std::cout);
        	bw_engine.write_dN_dt_table("tdndt.bin");
    	}
    	else{
            bw_engine.read_dN_dt_table("tdndt.bin");
    	}

    	if(!does_file_exist("tpair.bin")){
        	bw_engine.compute_cumulative_pair_table(&std::cout);
        	bw_engine.write_cumulative_pair_table("tpair.bin");
    	}
    	else{
        	bw_engine.read_cumulative_pair_table("tpair.bin");
    	}

	//Generates enough random numbers
	curandGenerator_t gen;
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, seed);
	double* d_rand;
	cudaMalloc(&d_rand, N*sizeof(double));
	curandGenerateUniformDouble(gen, d_rand, N);

	//Allocate on GPU an array of optical depths
	double* d_optical;
	cudaMalloc(&d_optical, N*sizeof(double));

	//Initialize the optical depths on GPU
	init_opt_depth<<<(N+255)/256, 256>>>(N, d_optical, d_rand);
	cudaDeviceSynchronize();

	//Copy back to the host & print for diag
	double* optical = new double[N];
	cudaMemcpy(optical, d_optical, N*sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "Test optical depths: " << std::endl;
	std::cout << optical[0] << " " << optical[1] << " " << optical[N/2] << " " << optical[N-2] << " " << optical[N-1] << std::endl;
	std::cout << "_________" << std::endl << std::endl;

	//Allocate space for momenta & fields & weigths.
	double* d_px;
	double* d_py;
	double* d_pz;
	double* d_ex;
	double* d_ey;
	double* d_ez;
	double* d_bx;
	double* d_by;
	double* d_bz;
	double* d_w;
	cudaMalloc(&d_px, N*sizeof(double));
	cudaMalloc(&d_py, N*sizeof(double));
	cudaMalloc(&d_pz, N*sizeof(double));
	cudaMalloc(&d_ex, N*sizeof(double));
	cudaMalloc(&d_ey, N*sizeof(double));
	cudaMalloc(&d_ez, N*sizeof(double));
	cudaMalloc(&d_bx, N*sizeof(double));
	cudaMalloc(&d_by, N*sizeof(double));
	cudaMalloc(&d_bz, N*sizeof(double));
	cudaMalloc(&d_w, N*sizeof(double));

	//Initialize momenta&fields randomly and print a test
	double* d_rand2;
	cudaMalloc(&d_rand2, N*sizeof(double)*9);
	curandGenerateUniformDouble(gen, d_rand2, N*9);
	init_mom_fields<<<(N+255)/256, 256>>>(N, d_px, d_py, d_pz, d_ex, d_ey, d_ez, d_bx, d_by, d_bz, d_w, d_rand2);
	cudaDeviceSynchronize();
	double px, py, pz, ex, ey, ez, bx, by, bz, w;
	cudaMemcpy(&px, &d_px[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&py, &d_py[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&pz, &d_pz[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&ex, &d_ex[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&ey, &d_ey[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&ez, &d_ez[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&bx, &d_bx[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&by, &d_by[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&bz, &d_bz[0], sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&w, &d_w[0], sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "Mom & fields & weight : " << std::endl;
	std::vector<double> vv{px,py,pz,ex,ey,ez,bx,by,bz,w};
	for (auto el: vv)
		std::cout << el << " ";
	std::cout << std::endl;
	std::cout << "_________" << std::endl << std::endl;

	//Export innards (in order to copy BW engin data to the GPU)
	auto innards = bw_engine.export_innards();

	//Copy TTfunc_table & bw_ctrl to GPU
	double* d_TTfunc_table_coords;
	double* d_TTfunc_table_data;
	cudaMalloc(&d_TTfunc_table_coords, sizeof(double)*innards.TTfunc_table_coords_how_many);
	cudaMalloc(&d_TTfunc_table_data, sizeof(double)*innards.TTfunc_table_coords_how_many);
	cudaMemcpy(d_TTfunc_table_coords, innards.TTfunc_table_coords_ptr,
		sizeof(double)*innards.TTfunc_table_coords_how_many, cudaMemcpyHostToDevice);
	cudaMemcpy(d_TTfunc_table_data, innards.TTfunc_table_data_ptr,
		sizeof(double)*innards.TTfunc_table_coords_how_many, cudaMemcpyHostToDevice);
	pxrmp::breit_wheeler_engine_ctrl<double>* d_bw_ctrl;
	cudaMalloc(&d_bw_ctrl, sizeof(pxrmp::breit_wheeler_engine_ctrl<double>));
	cudaMemcpy(d_bw_ctrl, &innards.bw_ctrl, sizeof(pxrmp::breit_wheeler_engine_ctrl<double>), cudaMemcpyHostToDevice);

	//Allocate space for has_event_happened and event_dt on GPU
	bool* d_has_event_happened;
	double* d_event_dt;
	cudaMalloc(&d_has_event_happened, sizeof(bool)*N);
	cudaMalloc(&d_event_dt, sizeof(double)*N);

	//Test internal_evolve_opt_depth_and_determine_event on GPU (multiple times!)
	for(int i = 0; i < repeat; i++){
		test_internal_evolve_opt_depth_and_determine_event<<<(N+255)/256, 256>>>
		(N, d_px, d_py, d_pz, d_ex, d_ey, d_ez, d_bx, d_by, d_bz, timestep, d_optical,
		d_has_event_happened, d_event_dt, innards.TTfunc_table_coords_how_many, d_TTfunc_table_coords,
		d_TTfunc_table_data, d_bw_ctrl);
		cudaDeviceSynchronize();
	}

	//Copy optical depths back to the host & print for diag
	double* optical2 = new double[N];
	cudaMemcpy(optical2, d_optical, N*sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "Test optical depths: " << std::endl;
	std::cout << optical[0] << " --> " << optical2[0] << std::endl;
	std::cout << optical[1] << " --> " << optical2[1] << std::endl;
	std::cout << optical[N/2] << " --> " << optical2[N/2] << std::endl;
	std::cout << optical[N-2] << " --> " << optical2[N-2] << std::endl;
	std::cout << optical[N-1] << " --> " << optical2[N-1] << std::endl;
	std::cout << "_________" << std::endl << std::endl;

	//For TEST PURPOSES, set all the elements of d_has_event_happened to TRUE
	set_all_events_true<<<(N+255)/256, 256>>>(N, d_has_event_happened);


	//Copy cum_distrib_table from the CPU to the GPU
	double* d_cum_distrib_table_coords_1;
	double* d_cum_distrib_table_coords_2;
	double* d_cum_distrib_table_data;
	cudaMalloc(&d_cum_distrib_table_coords_1, sizeof(double)*innards.cum_distrib_table_coords_1_how_many);
	cudaMalloc(&d_cum_distrib_table_coords_2, sizeof(double)*innards.cum_distrib_table_coords_2_how_many);
	cudaMalloc(&d_cum_distrib_table_data, sizeof(double)*innards.cum_distrib_table_coords_1_how_many*innards.cum_distrib_table_coords_2_how_many);
	cudaMemcpy(d_cum_distrib_table_coords_1, innards.cum_distrib_table_coords_1_ptr,
		sizeof(double)*innards.cum_distrib_table_coords_1_how_many,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cum_distrib_table_coords_2, innards.cum_distrib_table_coords_2_ptr,
		sizeof(double)*innards.cum_distrib_table_coords_2_how_many,cudaMemcpyHostToDevice);
	cudaMemcpy(d_cum_distrib_table_data, innards.cum_distrib_table_data_ptr,
		sizeof(double)*innards.cum_distrib_table_coords_1_how_many*innards.cum_distrib_table_coords_2_how_many,cudaMemcpyHostToDevice);

	//Generate enough random numbers on the GPU
	double* d_rand3;
	cudaMalloc(&d_rand3, sizeof(double)*N*sampling);
	curandGenerateUniformDouble(gen, d_rand3, N*sampling);

	//Allocate space for momenta & weigths of the generated pairs
	double* d_e_px;
	double* d_e_py;
	double* d_e_pz;
	double* d_p_px;
	double* d_p_py;
	double* d_p_pz;
	double* d_e_w;
	double* d_p_w;
	cudaMalloc(&d_e_px, sizeof(double)*N*sampling);
	cudaMalloc(&d_e_py, sizeof(double)*N*sampling);
	cudaMalloc(&d_e_pz, sizeof(double)*N*sampling);
	cudaMalloc(&d_p_px, sizeof(double)*N*sampling);
	cudaMalloc(&d_p_py, sizeof(double)*N*sampling);
	cudaMalloc(&d_p_pz, sizeof(double)*N*sampling);
	cudaMalloc(&d_e_w, sizeof(double)*N*sampling);
	cudaMalloc(&d_p_w, sizeof(double)*N*sampling);

	//Test internal_generate_breit_wheeler_pairs on the GPU
	test_internal_generate_breit_wheeler_pairs<<<(N+255)/256, 256>>>
	(N, d_has_event_happened,
	d_px, d_py, d_pz, d_ex, d_ey, d_ez, d_bx, d_by, d_bz,
	d_w, sampling,
	d_e_px, d_e_py, d_e_pz, d_p_px, d_p_py, d_p_pz, d_e_w, d_p_w,
	innards.cum_distrib_table_coords_1_how_many, d_cum_distrib_table_coords_1,
	innards.cum_distrib_table_coords_2_how_many, d_cum_distrib_table_coords_2, d_cum_distrib_table_data,
	d_bw_ctrl,
	d_rand3);


	//Copy some pair properties back to CPU for debug purposes
	double e_px, e_py, e_pz, e_w;
	double p_px, p_py, p_pz, p_w;
	cudaMemcpy(&e_px, d_e_px, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&e_py, d_e_py, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&e_pz, d_e_pz, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&e_w, d_e_w, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&p_px, d_p_px, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&p_py, d_p_py, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&p_pz, d_p_pz, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&p_w, d_p_w, sizeof(double), cudaMemcpyDeviceToHost);
	std::cout << "Test pairs: " << std::endl;
	std::cout << "e- : " << e_px << " " << e_py << " " << e_pz << " " << e_w << std::endl;
	std::cout << "e+ : " << p_px << " " << p_py << " " << p_pz << " " << p_w << std::endl;
	std::cout << "_________" << std::endl << std::endl;



	//Clean-up
	cudaFree(d_optical);
	cudaFree(d_rand);
	cudaFree(d_px);
	cudaFree(d_py);
	cudaFree(d_pz);
	cudaFree(d_ex);
	cudaFree(d_ey);
	cudaFree(d_ez);
	cudaFree(d_bx);
	cudaFree(d_by);
	cudaFree(d_bz);
	cudaFree(d_w);
	cudaFree(d_rand2);
	cudaFree(d_TTfunc_table_coords);
	cudaFree(d_TTfunc_table_data);
	cudaFree(d_bw_ctrl);
	cudaFree(d_has_event_happened);
	cudaFree(d_event_dt);
	cudaFree(d_rand3);
	cudaFree(d_e_px);
	cudaFree(d_e_py);
	cudaFree(d_e_pz);
	cudaFree(d_p_px);
	cudaFree(d_p_py);
	cudaFree(d_p_pz);
	cudaFree(d_e_w);
	cudaFree(d_p_w);
	delete[] optical;
	delete[] optical2;

	return 0;
}
