#include <iostream>
#include <array>

#include "quantum_sync_engine_tabulated_functions.hpp"

using namespace picsar::multi_physics::phys::quantum_sync;

int main()
{
	for (auto chi : std::array<double,5>{0.01, 0.1, 1.0, 10.0, 100.0}){
		for (auto xi : std::array<double,4>{0.01, 0.1, 0.5, 0.9}){
			std::cout << chi << " " << xi*chi << " " << compute_G_integrand<double>(chi, xi*chi) << std::endl;
		}
	}
		
	return 0;
}
