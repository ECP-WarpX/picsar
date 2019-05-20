//This header file contains usefuls definitions (macros, physical constants...)
//Should be included by all the source files of the QED library before

//The library allows to use either normalized or SI units.
//These compile-time directives protect against mixing
#if defined(PXRMP_USE_NORMALIZED_UNITS) && defined(PXRMP_USE_SI_UNITS)
  #error The library cannot be used mixing normalized and SI units!
#endif

#ifndef __PICSAR_MULTIPHYSICS_QED_COMMONS__
#define __PICSAR_MULTIPHYSICS_QED_COMMONS__

//###################### Decorator to compile some functions for CPU and GPU####

    //This flag should be set by the user
    //e.g. #define PXRMP_GPU __host__ __device__

    //If the user has not set the GPU flag
    #ifndef PXRMP_GPU
      //set it to the empty string
      #define PXRMP_GPU
    #else
      //otherwise, set this flag to inform the code that GPU code has to be
      //generated (this disables error messages...)
      #define PXRMP_WITH_GPU
    #endif

//##############################################################################


//############################################## Compiler specific macros ######

        //Restrict qualifier (compiler specific!)
        #ifdef _WIN32
            #define PXRMP_RESTRICT __restrict
        #else
            #define PXRMP_RESTRICT __restrict__
        #endif

        //Force inline pragmas (compiler specific!)
        #if defined(__CUDA_ARCH__)
          #define PXRMP_FORCE_INLINE __forceinline__
        #elif defined(__INTEL_COMPILER)
          #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__clang__)
          #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__GNUC__)
          #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined(__ibmxl__)
          #define PXRMP_FORCE_INLINE inline __attribute__((always_inline))
        #else
          #define PXRMP_FORCE_INLINE inline
        #endif


//##############################################################################

namespace picsar{
    namespace multi_physics{

//######################## Physical & mathematical constants (S.I.) ############

        //Double precision
        const double pi = 3.141592653589793;

        const double electron_mass = 9.10938356e-31;
        const double elementary_charge = 1.6021766208e-19;
        const double light_speed = 299792458.;
        const double reduced_plank = 1.054571800e-34;
        const double vacuum_permittivity =  8.854187817e-12;
        const double fine_structure =  0.0072973525664;

        const double classical_electron_radius =
        elementary_charge*elementary_charge /
        (4.0*pi*vacuum_permittivity*electron_mass*light_speed*light_speed);

        const double schwinger_field =
        electron_mass*electron_mass*(light_speed*light_speed*light_speed)/
        (elementary_charge*reduced_plank);

        const double pair_prod_rate_coeff =
        fine_structure * electron_mass * light_speed * light_speed /
        (reduced_plank);

        const double tau_e = classical_electron_radius/light_speed;

        const double quantum_synchrotron_rate_coeff =
        fine_structure*fine_structure/tau_e;

        const double si_gigameter = 1.0e9;
        const double si_megameter = 1.0e6;
        const double si_kilometer = 1.0e3;
        const double si_meter = 1.0;
        const double si_decimeter = 1.0e-1;
        const double si_centimeter = 1.0e-2;
        const double si_millimeter = 1.0e-3;
        const double si_micrometer = 1.0e-6;
        const double si_nanometer = 1.0e-9;
        const double si_picometer = 1.0e-12;
        const double si_femtometer = 1.0e-15;

//##############################################################################

//######################## Compile-time flags for SI or normalized units #######

//The library should never use directly physical constants, but rather the
//constants defined below, which take into account the normalization choice

    //Default is SI units
    #ifdef PXRMP_USE_NORMALIZED_UNITS
      #pragma message("Library will use normalized units")
      #define PXRMP_WITH_NORMALIZED_UNITS
    #else
      #pragma message("Library will use SI")
      #define PXRMP_WITH_SI_UNITS
    #endif

    // #ifdef PXRMP_USE_NORMALIZED_UNITS
    // #pragma message("Library will use normalized units")
    // #define PXRMP_WITH_SI_UNITS 0
    // #else
    // #pragma message("Library will use SI")
    // #define PXRMP_WITH_SI_UNITS 1
    // #endif

    #ifdef PXRMP_WITH_SI_UNITS
      const double __c = light_speed;
      const double __emass = electron_mass;
      const double __schwinger = schwinger_field;
      const double __pair_prod_coeff = pair_prod_rate_coeff*electron_mass*
        light_speed*light_speed;
      const double __quantum_synchrotron_rate_coeff =
        quantum_synchrotron_rate_coeff;
    #else
      const double __c = 1.0;
      const double __emass = 1.0;
      const double __schwinger = electron_mass*light_speed/(reduced_plank*2.*pi);
      const double __pair_prod_coeff = pair_prod_rate_coeff/(2.0*pi*light_speed);
      const double __quantum_synchrotron_rate_coeff =
        quantum_synchrotron_rate_coeff/(2.0*pi*light_speed);
    #endif

//##############################################################################


//######################## Default values for the parameters of BW engine ######

    //Minimum chi_phot to consider in all the functions
    //If chi is less than __bw_min_chi_phot
    //the function evolve_opt_depth_and_determine_event
    //provided by the BW engine will return immediately
    //without touching the optical depth and communicating that no event
    //has occurred.
    const double __breit_wheeler_min_chi_phot = 0.001;

    // dN/dt table:

    //__breit_wheeler_min_tdndt_chi_phot  is the inferior
    //limit of the total pair production rate lookup table.
    //If  __breit_wheeler_min_chi_phot < chi <__breit_wheeler_min_tdndt_chi_phot
    //BW process is taken into account, but the Erber approximation
    //rather than the lookup table is used.
    const double __breit_wheeler_min_tdndt_chi_phot = 0.1; //Min chi_phot
    const double __breit_wheeler_max_tdndt_chi_phot = 500.0;  //Max chi_phot
    const size_t __breit_wheeler_how_many_tdndt_chi_phot = 50;   //How many points
    // -------

    //Coefficients for che asymptotic behaviour of the "TT function"
    //using Erber approximation
    //Erber T. (1966), Reviews of Modern Physics,38, 626
    const double __erber_Tfunc_asynt_a = 0.16;
    const double __erber_Tfunc_asynt_b = 4./3.;

    // Pair production table:
    const double __breit_wheeler_min_tpair_chi_phot = 0.01; //Min chi_phot
    const double __breit_wheeler_max_tpair_chi_phot = 500.0; //Max chi_phot
    const size_t __breit_wheeler_how_many_tpair_chi_phot = 50; //How many points

    const size_t __breit_wheeler_chi_frac_tpair_how_many = 50;

    //Sets the limits for the semi-infinite integrals in the library
    const double __breit_wheeler_special_func_big_arg = 1.0e20;

//##############################################################################


//# Default values for the parameters of Quantum synchrotron engine ############

    //Minimum chi for particles to be considered by the engine

    const double __quantum_synchrotron_min_chi_part = 0.0001;

    // dN/dt table:

    // TO_WRITE
    const double __quantum_synchrotron_min_tdndt_chi_part = 0.0001;//Min chi_part
    const double __quantum_synchrotron_max_tdndt_chi_part = 200.0;//Max chi_part
    //How many points
    const size_t __quantum_synchrotron_how_many_tdndt_chi_part = 50;

    // Photon emission table:
    const double __quantum_synchrotron_min_tem_chi_part = 0.0001; //Min chi_part
    const double __quantum_synchrotron_max_tem_chi_part = 200.0; //Max chi_part
    const size_t __quantum_synchrotron_how_many_tem_chi_part = 50; //How many points

    const size_t __quantum_synchrotron_chi_frac_tem_how_many = 50;

    //Sets the limits for the semi-infinite integrals in the library
    const double __quantum_synchrotron_special_func_big_arg = 1.0e20;

// -------



//##############################################################################

    }
}

#endif// __PICSAR_MULTIPHYSICS_QED_COMMONS__
