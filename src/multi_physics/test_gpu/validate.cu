#include <iostream>
#include <vector>
#include <cstdio>
#include <cuda.h>
#include <curand.h>

#define PXRMP_GPU __host__ __device__
#define PXRMP_WITH_SI_UNITS
#include "../QED/src/breit_wheeler_engine.hpp"
#include "../QED/src/quantum_sync_engine.hpp"
#include "../QED/src/rng_wrapper.hpp"

//Alias for the picsar::multi_physics namespace
namespace pxrmp =  picsar::multi_physics;

//Seed of the random number generator
const size_t seed_bw = 83871734;
const size_t seed_qs = 93012221;
const size_t useless_seed = 22051988;

//A lot of particles!
const size_t N = 2000000;

//Physical constants
const double las_wavlngth = 1000.0 * pxrmp::si_nanometer;
const double me_c = pxrmp::electron_mass * pxrmp::light_speed;
const double eref = 2.0*pxrmp::pi*pxrmp::electron_mass*pxrmp::light_speed*pxrmp::light_speed/
            (las_wavlngth*pxrmp::elementary_charge);
const double bref = eref/pxrmp::light_speed;
const double one_femto = 1.0e-15;
//____________________________

//An empty class
//(will use RNG provided by cuda, not a wrapper of the
//RNG provided by the STL)
class dummy{};

//Lambda is not needed if SI units are used
const double default_lambda = 1.0;

//Helper function
bool does_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

void do_bw();

void do_qs();



int main()
{

    do_bw();
    do_qs();

}

__global__
void init_mom_fields(int n, double* px, double* py, double* pz,
	double*ex, double* ey, double* ez, double* bx, double* by, double* bz,
    double* ww, double* rand_nums,
    double min_px, double max_px,
    double min_py, double max_py,
    double min_pz, double max_pz,
    double min_ex, double max_ex,
    double min_ey, double max_ey,
    double min_ez, double max_ez,
    double min_bx, double max_bx,
    double min_by, double max_by,
    double min_bz, double max_bz)
{

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < n){
		px[i] = me_c*(min_px+(max_px-min_px)*rand_nums[0*n + i]);
		py[i] = me_c*(min_py+(max_py-min_py)*rand_nums[1*n + i]);
		pz[i] = me_c*(min_pz+(max_pz-min_pz)*rand_nums[2*n + i]);
		ex[i] = eref*(min_ex+(max_ex-min_ex)*rand_nums[3*n + i]);
		ey[i] = eref*(min_ey+(max_ey-min_ey)*rand_nums[4*n + i]);
		ez[i] = eref*(min_ez+(max_ez-min_ez)*rand_nums[5*n + i]);
		bx[i] = bref*(min_bx+(max_bx-min_bx)*rand_nums[6*n + i]);
		by[i] = bref*(min_by+(max_by-min_by)*rand_nums[7*n + i]);
		bz[i] = bref*(min_bz+(max_bz-min_bz)*rand_nums[8*n + i]);
		ww[i] = 1.0; // All weights are set to 1
	}
}



//*********************** BW ENGINE: test pair production rate ******************************
//GPU kernel to test internal_evolve_opt_depth_and_determine_event
__global__
void test_bw_prod(
	int n, double* px, double* py, double* pz, double* ex, double* ey, double* ez, double* bx, double* by, double* bz,
	size_t tab_how_many, double* coords, double* data, pxrmp::breit_wheeler_engine_ctrl<double>* bw_ctrl,
    double* chi, double* rate
)
{
	//Regenerate the lookuptable on GPU
	//This constructor does NOT allocate new memory: it manages existing pointers.
	pxrmp::lookup_1d<double> TTfunctab{tab_how_many, coords, data};

	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < n){
        double opt = 0.0;
        bool bdummy;
        double ddummy;
		pxrmp::breit_wheeler_engine<double, dummy>::
		internal_evolve_opt_depth_and_determine_event(
         		px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i],
         		one_femto, opt, bdummy, ddummy,
         		default_lambda, TTfunctab, *bw_ctrl);
        chi[i] = pxrmp::chi_photon<double>(px[i], py[i], pz[i],
            ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda);
        rate[i] = -opt;
	}
}
//********************************************************************************************



void do_bw()
{
    //Change default table parameters in order to speed up the calculations
    //pxrmp::breit_wheeler_engine_ctrl<double> bw_ctrl;
    //bw_ctrl.chi_phot_tdndt_how_many = 200;
    //bw_ctrl.chi_phot_tpair_how_many = 3;
    //bw_ctrl.chi_frac_tpair_how_many = 3;

    //Initialize the BW engine
    auto bw_engine =
        pxrmp::breit_wheeler_engine<double, pxrmp::stl_rng_wrapper>
        {std::move(pxrmp::stl_rng_wrapper{useless_seed}), default_lambda};//, bw_ctrl};

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

        //Initialize RNG on GPU
    	curandGenerator_t gen;
    	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    	curandSetPseudoRandomGeneratorSeed(gen, seed_bw);

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

        //Initialize momenta&fields
    	double* d_rand;
    	cudaMalloc(&d_rand, N*sizeof(double)*9);
    	curandGenerateUniformDouble(gen, d_rand, N*9);
        double pxmin = 0;
        double pymin = 0;
        double pzmin = -500000;
        double exmin = 500;
        double eymin = 0;
        double ezmin = 0;
        double bxmin = 0;
        double bymin = 0;
        double bzmin = 0;
        double pxmax = 0;
        double pymax = 0;
        double pzmax = 500000;
        double exmax = 500;
        double eymax = 0;
        double ezmax = 0;
        double bxmax = 0;
        double bymax = 0;
        double bzmax = 0;
    	init_mom_fields<<<(N+255)/256, 256>>>
        (N, d_px, d_py, d_pz, d_ex,
        d_ey, d_ez, d_bx, d_by, d_bz,
        d_w, d_rand,
        pxmin, pxmax,
        pymin, pymax,
        pzmin, pzmax,
        exmin, exmax,
        eymin, eymax,
        ezmin, ezmax,
        bxmin, bxmax,
        bymin, bymax,
        bzmin, bzmax);

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


        //test BW production rate & print on disk
        double* d_chi;
        double* d_rate;
        cudaMalloc(&d_chi, N*sizeof(double));
        cudaMalloc(&d_rate, N*sizeof(double));
        test_bw_prod<<<(N+255)/256, 256>>>(
        	N, d_px, d_py, d_pz, d_ex, d_ey, d_ez, d_bx, d_by, d_bz,
        	innards.TTfunc_table_coords_how_many, d_TTfunc_table_coords, d_TTfunc_table_data, d_bw_ctrl,
            d_chi, d_rate);
        double* chi = new double[N];
        double* rate = new double[N];
        cudaMemcpy(chi, d_chi, N*sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(rate, d_rate, N*sizeof(double), cudaMemcpyDeviceToHost);
        std::ofstream of{"bw_rate.dat"};
        for(size_t i = 0; i < N ; i++)
                of << chi[i] << " " << rate[i] << std::endl;

        of.close();



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
        cudaFree(d_rand);
        cudaFree(d_chi);
        cudaFree(d_rate);
        delete[] chi;
        delete[] rate;
}

void do_qs()
{

    //Change default table parameters in order to speed up the calculations
    //pxrmp::quantum_synchrotron_engine_ctrl<double> qs_ctrl;
    //qs_ctrl.chi_part_tdndt_how_many = 200;
    //qs_ctrl.chi_part_tem_how_many = 3;
    //qs_ctrl.chi_frac_tem_how_many = 3;

    //Initialize the QS engine
    auto qs_engine =
        pxrmp::quantum_synchrotron_engine<double, pxrmp::stl_rng_wrapper>
    	{std::move(pxrmp::stl_rng_wrapper{useless_seed}), default_lambda};//, qs_ctrl};

    if(!does_file_exist("qs_tdndt.bin")){
        qs_engine.compute_dN_dt_lookup_table(&std::cout);
        qs_engine.write_dN_dt_table("qs_tdndt.bin");
    }
    else{
        qs_engine.read_dN_dt_table("qs_tdndt.bin");
    }
    if(!does_file_exist("qs_photem.bin")){
        qs_engine.compute_cumulative_phot_em_table(&std::cout);
        qs_engine.write_cumulative_phot_em_table("qs_photem.bin");
    }
    else{
        qs_engine.read_cumulative_phot_em_table("qs_photem.bin");
    }

    //Initialize RNG on GPU
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, seed_qs);

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

    //Initialize momenta
    double* d_rand;
    cudaMalloc(&d_rand, N*sizeof(double)*9);
    curandGenerateUniformDouble(gen, d_rand, N*9);
    double pxmin = 0;
    double pymin = 0;
    double pzmin = -100000;
    double exmin = 500;
    double eymin = 0;
    double ezmin = 0;
    double bxmin = 0;
    double bymin = 0;
    double bzmin = 0;
    double pxmax = 0;
    double pymax = 0;
    double pzmax = 100000;
    double exmax = 500;
    double eymax = 0;
    double ezmax = 0;
    double bxmax = 0;
    double bymax = 0;
    double bzmax = 0;
    init_mom_fields<<<(N+255)/256, 256>>>
    (N, d_px, d_py, d_pz, d_ex,
    d_ey, d_ez, d_bx, d_by, d_bz,
    d_w, d_rand,
    pxmin, pxmax,
    pymin, pymax,
    pzmin, pzmax,
    exmin, exmax,
    eymin, eymax,
    ezmin, ezmax,
    bxmin, bxmax,
    bymin, bymax,
    bzmin, bzmax);

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
    cudaFree(d_rand);
}
