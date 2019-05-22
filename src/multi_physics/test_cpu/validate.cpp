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
const size_t N = 1000000;

//Physical constants
const double las_wavlngth = 1000.0 * pxrmp::si_nanometer;
const double me_c = pxrmp::electron_mass * pxrmp::light_speed;
const double eref = 2.0*pxrmp::pi*pxrmp::electron_mass*pxrmp::light_speed*pxrmp::light_speed/
            (las_wavlngth*pxrmp::elementary_charge);
const double bref = eref/pxrmp::light_speed;
const double one_femto = 1.0e-15;
const double me_c2 = pxrmp::electron_mass * pxrmp::light_speed * pxrmp::light_speed;
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

    //do_bw();
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

    auto innards = bw_engine.export_innards();
    auto c1 = innards.cum_distrib_table_coords_1_ptr;
    auto c2 = innards.cum_distrib_table_coords_2_ptr;
    auto dd = innards.cum_distrib_table_data_ptr;
    int cc = 0;
    std::vector<std::vector<double>> tab{
        innards.cum_distrib_table_coords_1_how_many,
        std::vector<double>(innards.cum_distrib_table_coords_2_how_many)};
    std::ofstream otab{"cum_prob_tab.dat"};
    for (int i = 0; i < innards.cum_distrib_table_coords_1_how_many; i++){
        for(int j = 0; j < innards.cum_distrib_table_coords_2_how_many; j++){
            tab[i][j] = exp(dd[cc]);
            otab << c1[i] << " " << c2[j] << " "  << exp(dd[cc++]) << std::endl;
        }
    }
    otab.close();
    std::ofstream otab2{"prob_tab.dat"};
    for (int i = 0; i < innards.cum_distrib_table_coords_1_how_many; i++){
        for(int j = 1; j < innards.cum_distrib_table_coords_2_how_many; j++){
            otab2 << exp(c1[i]) << " " << c2[j] << " "  << (tab[i][j]-tab[i][j-1]) << std::endl;
        }
    }
    otab2.close();

    //Initialize momenta&fields
    pxmin = 1500;
    pymin = 0;
    pzmin = 0;
    exmin = 0;
    eymin = 0;
    ezmin = 0;
    bxmin = 0;
    bymin = 0;
    bzmin = 270;
    pxmax = 1500;
    pymax = 0;
    pzmax = 0;
    exmax = 0;
    eymax = 0;
    ezmax = 0;
    bxmax = 0;
    bymax = 0;
    bzmax = 270;
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

    //test BW pair properties & print on disk
    std::vector<double> gamma_ele(N);
    std::vector<double> gamma_pos(N);
    for(size_t i = 0; i < N; i++){
        auto all = bw_engine.generate_breit_wheeler_pairs(
            px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i], w[i], 1);
        auto e_mom = all[0][0].first;
        auto p_mom = all[1][0].first;

        double norm_pe2 = pxrmp::norm2(e_mom);
        double norm_pp2 = pxrmp::norm2(p_mom);
        gamma_ele[i] = sqrt(1.0 + norm_pe2/me_c/me_c);
        gamma_pos[i] = sqrt(1.0 + norm_pp2/me_c/me_c);
    }

    std::ofstream offf{"smicheck.dat"};
    for(size_t i = 0; i < N; i++)
        offf << gamma_ele[i] << " " << gamma_pos[i] << std::endl;
    offf.close();


}

void do_qs()
{
    //Initialize the BW engine
    auto qs_engine =
        pxrmp::quantum_synchrotron_engine<double, pxrmp::stl_rng_wrapper>
        {std::move(pxrmp::stl_rng_wrapper{seed_qs}), default_lambda};//, bw_ctrl};

        //Initialize the lookup tables
        //Generates tables if they do not exist
    if(!does_file_exist("em_tdndt.bin")){
        qs_engine.compute_dN_dt_lookup_table(&std::cout);
        qs_engine.write_dN_dt_table("em_tdndt.bin");
    }
    else{
        qs_engine.read_dN_dt_table("em_tdndt.bin");
    }
    if(!does_file_exist("em_tphot.bin")){
        qs_engine.compute_cumulative_phot_em_table(&std::cout);
        qs_engine.write_cumulative_phot_em_table("em_tphot.bin");
    }
    else{
        qs_engine.read_cumulative_phot_em_table("em_tphot.bin");
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
    double pxmin = -50000;
    double pymin = 0;
    double pzmin = 0;
    double exmin = 0;
    double eymin = 0;
    double ezmin = 0;
    double bxmin = 0;
    double bymin = 0;
    double bzmin = 1000;
    double pxmax = 50000;
    double pymax = 0;
    double pzmax = 0;
    double exmax = 0;
    double eymax = 0;
    double ezmax = 0;
    double bxmax = 0;
    double bymax = 0;
    double bzmax = 1000;
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

        //test QS production rate & print on disk
    std::vector<double> chi(N);
    std::vector<double> rate(N);
    for(size_t i = 0; i < N; i++){
            double opt = 0.0;
        qs_engine.evolve_opt_depth_and_determine_event(
                px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i],
                one_femto, opt);
            chi[i] = pxrmp::chi_lepton<double>(px[i], py[i], pz[i],
                    ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda);
            rate[i] = -opt;
    }
    std::ofstream of{"qs_rate.dat"};
        for(size_t i = 0; i < N ; i++)
                of << chi[i] << " " << rate[i] << std::endl;
    of.close();


    //Initialize momenta&fields
 pxmin = 1000;
 pymin = 0;
 pzmin = 0;
 exmin = 0;
 eymin = 0;
 ezmin = 0;
 bxmin = 0;
 bymin = 0;
 bzmin = 1000;
 pxmax = 1000;
 pymax = 0;
 pzmax = 0;
 exmax = 0;
 eymax = 0;
 ezmax = 0;
 bxmax = 0;
 bymax = 0;
 bzmax = 1000;
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


    //test QS photon properties & print on disk
    std::vector<double> gamma_phot(N);
    for(size_t i = 0; i < N; i++){
        chi[i] = pxrmp::chi_lepton<double>(px[i], py[i], pz[i],
                ex[i], ey[i], ez[i], bx[i], by[i], bz[i], default_lambda);
        auto phot = qs_engine.generate_photons_and_update_momentum(
    			px[i], py[i], pz[i], ex[i], ey[i], ez[i], bx[i], by[i], bz[i], w[i], 1);

        auto p_phot = phot[0].first;

        gamma_phot[i] = pxrmp::norm(p_phot)/me_c;
    }




    std::ofstream of_photen{"qs_photen.dat"};
    for(size_t i = 0; i < N ; i++)
            of_photen << gamma_phot[i] << std::endl;
    of_photen.close();


    auto innards = qs_engine.export_innards();
    auto cc = innards.KKfunc_table_coords_ptr;
    auto dd = innards.KKfunc_table_data_ptr;
    std::ofstream emtab1{"em_dndt_tab.dat"};
    for (int i = 0; i < innards.KKfunc_table_coords_how_many; i++){
            emtab1 << cc[i] << " " << exp(dd[i]) <<  std::endl;
    }
    emtab1.close();

    auto cc1 = innards.cum_distrib_table_coords_1_ptr;
    auto cc2 = innards.cum_distrib_table_coords_2_ptr;
    auto dd2 = innards.cum_distrib_table_data_ptr;

    std::ofstream emtab2{"em_photem.dat"};
    int count = 0;
    for (int i = 0; i < innards.cum_distrib_table_coords_1_how_many; i++){
        for (int j = 0; j < innards.cum_distrib_table_coords_2_how_many; j++){
            emtab2 << exp(cc1[i]) << " " << cc2[j] <<  " " << exp(dd2[count++]) <<  std::endl;
        }
    }
    emtab2.close();




}
