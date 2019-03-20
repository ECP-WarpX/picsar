#include <iostream>
#include <fstream>

//Import BW engine with SI units
#define PXRMP_USE_SI_UNITS
#include "breit_wheeler_engine.hpp"

//Defines an alias for the nested namespace
namespace pxrmp =  picsar::multi_physics;

int main(){
	pxrmp::stl_rng_wrapper wrap{22051988}; //Creates a wrapper around the STL RNG

	double laser_lambda = 1.0; //For SI units this paramter is ignored

	//This is a struct which contains parameters relevant for BW engine
	pxrmp::breit_wheeler_engine_ctrl<double> bw_ctrl;
	//sets the minimum chi for which calculations are actually done
	bw_ctrl.chi_phot_min = 0.01;

	//Creates the BW engine
	auto bw_engine =  pxrmp::breit_wheeler_engine<double, pxrmp::stl_rng_wrapper>{
		std::move(wrap),
		laser_lambda, //Optional argument (default value is 1)
		bw_ctrl //Optional argument (with reasonable default values)
		};

	//***SI units
	double lambda = 800.0*pxrmp::si_nanometer;
	const double me_c = pxrmp::electron_mass * pxrmp::light_speed;
	double eref = 2.0*pxrmp::pi*pxrmp::electron_mass*pxrmp::light_speed*
		pxrmp::light_speed/(lambda*pxrmp::elementary_charge);
	double bref = eref/pxrmp::light_speed;
	double dtref = lambda/(2.0*pxrmp::pi*pxrmp::light_speed);
	//double rateref = 1.0/dtref;
	//***

	//***Initial properties of a particle
      	double px =  149.825*me_c;
    	double py =  933.115*me_c;
    	double pz =  -538.195*me_c;
    	double ex =  931.686*eref;
    	double ey =  -861.074*eref;
    	double ez =  944.652*eref;
    	double bx =  531.406*bref;
    	double by =  670.933*bref;
    	double bz =  660.057*bref;
	//***

    //Timestep
	double dt = 0.01*dtref;

	//Initial optical depth
	double opt_depth = 0.001;

	bw_engine.compute_dN_dt_lookup_table(&std::cout);

	//Decrease optical depth & determine if an event has occured
	bool has_event_happend;
    	double dt_prod;
   	std::tie(has_event_happend, dt_prod) =
        bw_engine.evolve_opt_depth_and_determine_event
        (px, py, pz, ex, ey, ez, bx, by, bz, dt, opt_depth);

	if(has_event_happend){
		std::cout << "An event occured after: "  << dt_prod <<
			"s ( " << 100.0* dt_prod/dt << " %% of dt)" << std::endl;
	}


	std::cout << "Table vs direct calculation " << std::endl;
	std::cout << "chi	interp	calc" << std::endl;

	double Efield = 500.0 * eref;

	double mom = 0.01*me_c;
	double mul = 1.2;

	while(mom < 1e7 * me_c){
		double chi = pxrmp::chi_photon(mom,0.,0.,0.,Efield,0.,0.,0.,0.);
		double ephot = mom*pxrmp::light_speed;
		std::cout << chi << "   ";
		std::cout <<  bw_engine.interp_dN_dt(ephot, chi) << "   ";
		std::cout <<  bw_engine.compute_dN_dt(ephot, chi) << "\n";
		mom*=mul;
	}

	std::cout << "Saving cumulative pair production on disk..." << std::endl;

	std::ofstream of;
	of.open("cum_pair_prod.dat");

	double chi = 0.1;
	double mulc = 1.1;
	while(chi < 1.0e3){
		std::cout << "*" << chi << std::endl;
		for(int t = 1; t < 100; t++){
			double frac = t/100.0;
			of << chi << " " << frac << " ";
			of << bw_engine.compute_cumulative_pair(chi, chi*frac) << "\n";
		}
		of <<  "\n";	
		chi *=mulc;
	}



	of.close();


	return 0;
}
