#include <iostream>
#include <vector>
#include <random>

#include "../QED/src/breit_wheeler_engine.hpp"
#include "../QED/src/quantum_sync_engine.hpp"
#include "../QED/src/rng_wrapper.hpp"

//Alias for the picsar::multi_physics namespace
namespace pxrmp =  picsar::multi_physics;

//Seed of the random number generator
const size_t seed_bw = 83871734;
const size_t seed_qs = 93012221;
const size_t seed = 22051988;

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

void init_mom_fields(std::vector<double>& px, std::vector<double>& py, std::vector<double>& pz,
	std::vector<double>& ex, std::vector<double>& ey, std::vector<double>& ez, std::vector<double>& bx, std::vector<double>& by, std::vector<double>& bz,
    std::vector<double>& ww, std::default_random_engine rng,
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

  	std::uniform_real_distribution<double> distribution(0.0,1.0);
	for (size_t i = 0; i < px.size(); i++){
		px[i] = me_c*(min_px+(max_px-min_px)*distribution(rng));
		py[i] = me_c*(min_py+(max_py-min_py)*distribution(rng));
		pz[i] = me_c*(min_pz+(max_pz-min_pz)*distribution(rng));
		ex[i] = eref*(min_ex+(max_ex-min_ex)*distribution(rng));
		ey[i] = eref*(min_ey+(max_ey-min_ey)*distribution(rng));
		ez[i] = eref*(min_ez+(max_ez-min_ez)*distribution(rng));
		bx[i] = bref*(min_bx+(max_bx-min_bx)*distribution(rng));
		by[i] = bref*(min_by+(max_by-min_by)*distribution(rng));
		bz[i] = bref*(min_bz+(max_bz-min_bz)*distribution(rng));
		ww[i] = 1.0; // All weights are set to 1
	}
}

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
        {std::move(pxrmp::stl_rng_wrapper{seed_bw}), default_lambda};//, bw_ctrl};

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

        //Allocate space for momenta & fields & weigths.
        std::vector<double> px(N);
	std::vector<double> py(N);
 	std::vector<double> pz(N);
 	std::vector<double> ex(N);
 	std::vector<double> ey(N);
	std::vector<double> ez(N);
	std::vector<double> bx(N);
 	std::vector<double> by(N);
	std::vector<double> bz(N);
	std::vector<double> w(N);


        //Initialize momenta&fields
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
    	init_mom_fields
        (px, py, pz, ex,
        ey, ez, bx, by, bz,
        w, std::default_random_engine(seed),
        pxmin, pxmax,
        pymin, pymax,
        pzmin, pzmax,
        exmin, exmax,
        eymin, eymax,
        ezmin, ezmax,
        bxmin, bxmax,
        bymin, bymax,
        bzmin, bzmax);

	//test BW production rate & print on disk
	std::vector<double> chi(N);
	std::vector<double> rate(N);
	for(size_t i = 0; i < N; i++){
			double opt = 0.0;
		bw_engine.evolve_opt_depth_and_determine_event(
         		px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i],
         		one_femto, opt);
        	chi[i] = pxrmp::chi_photon<double>(px[i], py[i], pz[i],
            		ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda);
        	rate[i] = -opt;
	}
	std::ofstream of{"bw_rate.dat"};
        for(size_t i = 0; i < N ; i++)
                of << chi[i] << " " << rate[i] << std::endl;

        of.close();

	//test BW pair properties & print on disk
	std::vector<double> chi_ele_frac(N);
	std::vector<double> chi_pos_frac(N);
	for(size_t i = 0; i < N; i++){
		auto all = bw_engine.generate_breit_wheeler_pairs(
			px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i], w[i], 1);
		auto e_mom = all[0][0].first;
		auto p_mom = all[1][0].first;

		chi[i] = pxrmp::chi_photon<double>(px[i], py[i], pz[i],
                ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda);
            chi_ele_frac[i] =  pxrmp::chi_lepton<double>(e_mom[0], e_mom[1], e_mom[2],
                ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda)/chi[i];
            chi_pos_frac[i] =  pxrmp::chi_lepton<double>(p_mom[0], p_mom[1], p_mom[2],
                ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda)/chi[i];
	}

	std::ofstream of2{"bw_pairs.dat"};
        for(size_t i = 0; i < N ; i++)
                    of2 << chi[i] << " " << chi_ele_frac[i] << " " << chi_pos_frac[i] << std::endl;

        of2.close();



}

void do_qs()
{}
