! ________________________________________________________________________________________
!
! *** Copyright Notice ***
!
! "Particle In Cell Scalable Application Resource (PICSAR) v2", Copyright (c)
! 2016, The Regents of the University of California, through Lawrence Berkeley
! National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
!
! NOTICE.
! This Software was developed under funding from the U.S. Department of Energy
! and the U.S. Government consequently retains certain rights. As such, the U.S.
! Government has been granted for itself and others acting on its behalf a
! paid-up, nonexclusive, irrevocable, worldwide license in the Software to
! reproduce, distribute copies to the public, prepare derivative works, and
! perform publicly and display publicly, and to permit other to do so.
!
! NCI_GODFREY_CORR.F90
!
! Purpose:
! This file contains a module for the coefficients and the routines for corrector of 
! Numerical Cherenkov Radiation in the Yee/Cole-Karkkainen solver. It assumes the plasma 
! drifts along z with gamma >> 1
! See article Journal of Computational Physics 267 (2014) 1–6
! Suppressing the numerical Cherenkov instability in FDTD PIC codes
! BB Godfrey, JL Vay
! Authors:
! Henri Vincenti
! Maxence Thevenet
!
! Date:
! Creation 2018
! ________________________________________________________________________________________


! ________________________________________________________________________________________
!> @brief
!> Module containing tabulated stencil values for the corrector as well as 2 subroutines:
!> init_godfrey_filter_coeffs to get the stencil value from the table
!> apply_filter_z_2d to apply the filter at each timestep
!> It assumes the plasma drifts along z with gamma >> 1
!>
!> @author
!> Maxence Thevenet
!> Henri Vincenti
!> 
!> @date
!> Creation 2018
! ________________________________________________________________________________________

MODULE godfrey_filter_coeffs
  use iso_c_binding
  USE PICSAR_precision
  USE constants
IMPLICIT NONE
  REAL(num), DIMENSION(4,100), PARAMETER :: coeff_ex_galerkin = reshape([&
    -2.47536,2.04288,-0.598163,0.0314711, & 
    -2.47536,2.04288,-0.598163,0.0314711, & 
    -2.47545,2.04309,-0.598307,0.0315029, & 
    -2.4756,2.04342,-0.598549,0.0315558, & 
    -2.47581,2.0439,-0.598886,0.0316298, & 
    -2.47608,2.0445,-0.59932,0.031725, & 
    -2.47641,2.04525,-0.59985,0.0318412, & 
    -2.4768,2.04612,-0.600477,0.0319785, & 
    -2.47725,2.04714,-0.6012,0.0321367, & 
    -2.47776,2.04829,-0.602019,0.0323158, & 
    -2.47833,2.04957,-0.602934,0.0325158, & 
    -2.47896,2.05099,-0.603944,0.0327364, & 
    -2.47965,2.05254,-0.605051,0.0329777, & 
    -2.4804,2.05423,-0.606253,0.0332396, & 
    -2.48121,2.05606,-0.60755,0.0335218, & 
    -2.48208,2.05802,-0.608942,0.0338243, & 
    -2.48301,2.06012,-0.610429,0.0341469, & 
    -2.48401,2.06235,-0.61201,0.0344895, & 
    -2.48506,2.06471,-0.613685,0.0348519, & 
    -2.48618,2.06721,-0.615453,0.0352339, & 
    -2.48735,2.06984,-0.617314,0.0356353, & 
    -2.48859,2.07261,-0.619268,0.0360559, & 
    -2.48988,2.0755,-0.621312,0.0364954, & 
    -2.49123,2.07853,-0.623447,0.0369536, & 
    -2.49265,2.08169,-0.625672,0.0374302, & 
    -2.49412,2.08498,-0.627986,0.0379248, & 
    -2.49565,2.0884,-0.630386,0.0384372, & 
    -2.49724,2.09194,-0.632873,0.0389669, & 
    -2.49888,2.09561,-0.635443,0.0395135, & 
    -2.50058,2.09939,-0.638096,0.0400766, & 
    -2.50234,2.1033,-0.640829,0.0406557, & 
    -2.50415,2.10732,-0.64364,0.0412502, & 
    -2.50601,2.11145,-0.646526,0.0418594, & 
    -2.50791,2.1157,-0.649485,0.0424828, & 
    -2.50987,2.12004,-0.652512,0.0431196, & 
    -2.51187,2.12448,-0.655604,0.0437688, & 
    -2.51392,2.12901,-0.658756,0.0444297, & 
    -2.516,2.13363,-0.661964,0.0451011, & 
    -2.51812,2.13832,-0.665221,0.0457818, & 
    -2.52027,2.14308,-0.668521,0.0464705, & 
    -2.52244,2.14789,-0.671856,0.0471658, & 
    -2.52464,2.15274,-0.675218,0.0478658, & 
    -2.52684,2.15762,-0.678596,0.0485687, & 
    -2.52906,2.16251,-0.68198,0.0492723, & 
    -2.53126,2.16738,-0.685355,0.049974, & 
    -2.53345,2.17222,-0.688706,0.0506708, & 
    -2.53561,2.177,-0.692015,0.0513594, & 
    -2.53773,2.18168,-0.69526,0.0520359, & 
    -2.53978,2.18623,-0.698416,0.0526955, & 
    -2.54175,2.19059,-0.701452,0.053333, & 
    -2.5436,2.19471,-0.704331,0.0539417, & 
    -2.54531,2.19852,-0.70701,0.0545141, & 
    -2.54683,2.20193,-0.709433,0.0550409, & 
    -2.5481,2.20483,-0.711533,0.0555106, & 
    -2.54906,2.20709,-0.713224,0.0559094, & 
    -2.54963,2.20852,-0.714397,0.0562198, & 
    -2.54968,2.20888,-0.714907,0.0564196, & 
    -2.54905,2.20785,-0.714562,0.0564797, & 
    -2.54751,2.20496,-0.713094,0.0563618, & 
    -2.54472,2.19955,-0.710118,0.0560124, & 
    -2.54014,2.19058,-0.705048,0.0553544, & 
    -2.53286,2.1763,-0.69693,0.0542684, & 
    -2.52115,2.15344,-0.684027,0.05255, & 
    -2.50098,2.11466,-0.66255,0.0497817, & 
    -2.45797,2.03459,-0.620099,0.0446889, & 
    -2.28371,1.72254,-0.465905,0.0283268, & 
    -2.4885,2.04899,-0.599292,0.0390466, & 
    -2.1433,1.36735,-0.220924,-0.00215633, & 
    -2.4943,2.07019,-0.610552,0.035166, & 
    -2.84529,2.77303,-1.00018,0.0724884, & 
    -2.72242,2.51888,-0.847226,0.0509964, & 
    -2.65633,2.3744,-0.750392,0.0326366, & 
    -2.59601,2.23412,-0.646421,0.00868027, & 
    -2.51477,2.0369,-0.491066,-0.0306397, & 
    -2.35935,1.65155,-0.178971,-0.112713, & 
    -1.84315,0.361693,0.876104,-0.393844, & 
    -2.65422,2.39262,-0.789663,0.0516265, & 
    -3.46529,4.42354,-2.45543,0.497097, & 
    -3.15747,3.65311,-1.824,0.328432, & 
    -3.04694,3.37613,-1.59668,0.267631, & 
    -2.99205,3.23814,-1.48302,0.237103, & 
    -2.96075,3.15894,-1.41733,0.219317, & 
    -2.94172,3.11028,-1.37649,0.20811, & 
    -2.92994,3.07962,-1.35025,0.200755, & 
    -2.92283,3.06054,-1.33338,0.195859, & 
    -2.91894,3.04938,-1.3229,0.192637, & 
    -2.91736,3.04394,-1.31702,0.190612, & 
    -2.91753,3.04278,-1.31456,0.189477, & 
    -2.91905,3.04494,-1.31475,0.189026, & 
    -2.92165,3.04973,-1.31705,0.189117, & 
    -2.92512,3.05667,-1.32105,0.189646, & 
    -2.92933,3.06539,-1.32646,0.190538, & 
    -2.93416,3.07562,-1.33308,0.191735, & 
    -2.93952,3.08715,-1.34072,0.193194, & 
    -2.94535,3.09982,-1.34925,0.194881, & 
    -2.95159,3.11349,-1.35858,0.196769, & 
    -2.9582,3.12805,-1.36861,0.198838, & 
    -2.96514,3.14342,-1.37929,0.201068, & 
    -2.97239,3.15953,-1.39055,0.203448, & 
    -2.97991,3.17632,-1.40234,0.205964, & 
    -2.98769,3.19374,-1.41463,0.208607 ], [4,100])
  REAL(num), DIMENSION(4,100), PARAMETER :: coeff_by_galerkin = reshape([&
    -2.80862,2.80104,-1.14615,0.154077, & 
    -2.80862,2.80104,-1.14615,0.154077, & 
    -2.80851,2.80078,-1.14595,0.154027, & 
    -2.80832,2.80034,-1.14561,0.153945, & 
    -2.80807,2.79973,-1.14514,0.153829, & 
    -2.80774,2.79894,-1.14454,0.15368, & 
    -2.80733,2.79798,-1.1438,0.153498, & 
    -2.80685,2.79685,-1.14292,0.153284, & 
    -2.8063,2.79554,-1.14192,0.153036, & 
    -2.80568,2.79405,-1.14077,0.152756, & 
    -2.80498,2.79239,-1.1395,0.152443, & 
    -2.80421,2.79056,-1.13809,0.152098, & 
    -2.80337,2.78856,-1.13656,0.151721, & 
    -2.80246,2.78638,-1.13488,0.151312, & 
    -2.80147,2.78404,-1.13308,0.150871, & 
    -2.80041,2.78152,-1.13115,0.150397, & 
    -2.79927,2.77882,-1.12908,0.149893, & 
    -2.79807,2.77596,-1.12689,0.149356, & 
    -2.79679,2.77292,-1.12456,0.148789, & 
    -2.79543,2.76972,-1.12211,0.14819, & 
    -2.79401,2.76634,-1.11953,0.14756, & 
    -2.79251,2.76279,-1.11681,0.1469, & 
    -2.79094,2.75907,-1.11397,0.146208, & 
    -2.78929,2.75517,-1.111,0.145486, & 
    -2.78757,2.7511,-1.10789,0.144733, & 
    -2.78578,2.74686,-1.10466,0.14395, & 
    -2.78391,2.74245,-1.1013,0.143137, & 
    -2.78196,2.73786,-1.09781,0.142293, & 
    -2.77994,2.73309,-1.09419,0.141419, & 
    -2.77784,2.72814,-1.09043,0.140514, & 
    -2.77566,2.72301,-1.08654,0.139578, & 
    -2.7734,2.7177,-1.08252,0.138612, & 
    -2.77106,2.7122,-1.07836,0.137614, & 
    -2.76864,2.70651,-1.07406,0.136586, & 
    -2.76613,2.70062,-1.06962,0.135525, & 
    -2.76353,2.69453,-1.06503,0.134432, & 
    -2.76084,2.68824,-1.0603,0.133307, & 
    -2.75806,2.68173,-1.05541,0.132148, & 
    -2.75518,2.675,-1.05037,0.130954, & 
    -2.75219,2.66804,-1.04516,0.129725, & 
    -2.7491,2.66084,-1.03978,0.12846, & 
    -2.7459,2.65339,-1.03423,0.127156, & 
    -2.74257,2.64566,-1.02848,0.125813, & 
    -2.73912,2.63765,-1.02254,0.124428, & 
    -2.73552,2.62934,-1.01638,0.122999, & 
    -2.73178,2.62069,-1.01,0.121523, & 
    -2.72787,2.61169,-1.00337,0.119996, & 
    -2.72379,2.6023,-0.996479,0.118417, & 
    -2.71951,2.59248,-0.989294,0.116778, & 
    -2.71501,2.58218,-0.981786,0.115076, & 
    -2.71026,2.57135,-0.97392,0.113303, & 
    -2.70524,2.55991,-0.965651,0.111453, & 
    -2.69989,2.54778,-0.956922,0.109514, & 
    -2.69416,2.53484,-0.947666,0.107476, & 
    -2.68799,2.52096,-0.937795,0.105324, & 
    -2.68129,2.50596,-0.927197,0.103039, & 
    -2.67394,2.48959,-0.915724,0.100597, & 
    -2.66578,2.47153,-0.903179,0.097968, & 
    -2.65657,2.4513,-0.889283,0.0951084, & 
    -2.64598,2.42824,-0.873638,0.0919592, & 
    -2.63347,2.40127,-0.855632,0.0884325, & 
    -2.61813,2.36864,-0.834261,0.0843898, & 
    -2.59821,2.32701,-0.807691,0.0795876, & 
    -2.56971,2.26887,-0.77188,0.0735132, & 
    -2.51823,2.16823,-0.713448,0.0645399, & 
    -2.33537,1.8294,-0.533852,0.0409941, & 
    -2.53143,2.14818,-0.670502,0.053982, & 
    -2.17737,1.43641,-0.259095,0.00101255, & 
    -2.51929,2.12931,-0.654743,0.0452381, & 
    -2.86122,2.82221,-1.05039,0.0894636, & 
    -2.72908,2.54506,-0.87834,0.0626188, & 
    -2.6536,2.37954,-0.7665,0.0409117, & 
    -2.58374,2.21923,-0.649738,0.0146791, & 
    -2.49284,2.00346,-0.48457,-0.0255348, & 
    -2.32762,1.60337,-0.1698,-0.105287, & 
    -1.80149,0.316787,0.855414,-0.369652, & 
    -2.60242,2.28418,-0.721378,0.040091, & 
    -3.40335,4.25157,-2.29817,0.449834, & 
    -3.0852,3.47341,-1.67791,0.28982, & 
    -2.9642,3.17856,-1.44399,0.229852, & 
    -2.89872,3.01966,-1.31861,0.197945, & 
    -2.85668,2.91811,-1.23894,0.17783, & 
    -2.82679,2.84621,-1.18287,0.163785, & 
    -2.80401,2.79167,-1.14058,0.153278, & 
    -2.78577,2.74819,-1.10706,0.145015, & 
    -2.77061,2.7122,-1.07946,0.138267, & 
    -2.75764,2.68152,-1.05606,0.132589, & 
    -2.74627,2.65475,-1.03575,0.127695, & 
    -2.73612,2.63093,-1.01777,0.123395, & 
    -2.72692,2.6094,-1.00159,0.119553, & 
    -2.71846,2.58968,-0.986841,0.116074, & 
    -2.71061,2.57142,-0.973239,0.112887, & 
    -2.70323,2.55434,-0.960573,0.109937, & 
    -2.69626,2.53824,-0.948678,0.107185, & 
    -2.68962,2.52294,-0.937429,0.104598, & 
    -2.68327,2.50833,-0.926722,0.102151, & 
    -2.67714,2.4943,-0.916477,0.0998223, & 
    -2.67122,2.48076,-0.906627,0.0975966, & 
    -2.66546,2.46764,-0.897118,0.0954599, & 
    -2.65985,2.45489,-0.887903,0.0934011, & 
    -2.65437,2.44244,-0.878945,0.0914107 ], [4,100])

CONTAINS

! ________________________________________________________________________________________
!> @brief
!> Subroutine to initialize stencils for the Godfrey FDTD NCI corrector, assuming the
!> plasma is drifting at relativistic speed in the z direction. These coefficients
!> only depend on vz*dt/dz, approximated by c*dt/dz
!> 
!> @param[inout] stencilz_ex     stencil along z only for Ex, Ey, Bz
!> @param[inout] stencilz_by     stencil along z only for Ez, Bx, By
!> @param[in] nstencilz          stencil order. 5 for this NCI corrector
!> @param[in] cdtodz             c*dt/dz. determines the lines read in the table
!> @param[in] l_lower_order_in_v integer (0 or 1): current deposition performed with lower
!>                                   order shape factor in the direction of each deposited 
!>                                   component of J (i.e., Galerkin scheme). Other
!>                                   schemes are not implemented for the moment, so
!>                                   l_lower_order_in_v should be 1.
!>
!> @author
!> Maxence Thevenet
!> @date
!> Creation 2018
! ________________________________________________________________________________________
SUBROUTINE init_godfrey_filter_coeffs(stencilz_ex, stencilz_by, nstencilz, cdtodz, l_lower_order_in_v) &
  bind(c, name='init_godfrey_filter_coeffs')
  INTEGER, value, INTENT(IN) :: nstencilz
  REAL(num), value, INTENT(IN) :: cdtodz
  INTEGER, value, INTENT(IN) :: l_lower_order_in_v
  REAL(num), INTENT(INOUT) :: stencilz_ex(0:nstencilz-1), stencilz_by(0:nstencilz-1)
  REAL(num), DIMENSION(0:3) :: prestencil_ex, prestencil_by
  INTEGER :: index, i, size_coeff_table
  REAL(num) :: weight_right
  
  IF (nstencilz .NE. 5) THEN
    WRITE(*,*) "ERROR: NCI corrector must be initialized with nstencilz=5"
    STOP  
  ENDIF
  IF (l_lower_order_in_v .NE. 0) THEN
    size_coeff_table = SIZE(coeff_ex_galerkin, 2)
    index = INT(size_coeff_table*cdtodz)+1
    IF (index<1) THEN
      index = 1
    ENDIF
    IF (index>size_coeff_table-1) THEN
      index = size_coeff_table-1
    ENDIF
    weight_right = cdtodz - (index-1)/size_coeff_table
    DO i=0, 3
      prestencil_ex(i) = (1-weight_right)*coeff_ex_galerkin(i+1,index) + weight_right*coeff_ex_galerkin(i+1,index+1)
      prestencil_by(i) = (1-weight_right)*coeff_by_galerkin(i+1,index) + weight_right*coeff_by_galerkin(i+1,index+1)
    END DO
    stencilz_ex(0) =  (256+128*prestencil_ex(0)+96*prestencil_ex(1)+80*prestencil_ex(2)+70*prestencil_ex(3))/256
    stencilz_ex(1) = -(     64*prestencil_ex(0)+64*prestencil_ex(1)+60*prestencil_ex(2)+56*prestencil_ex(3))/256
    stencilz_ex(2) =  (                         16*prestencil_ex(1)+24*prestencil_ex(2)+28*prestencil_ex(3))/256
    stencilz_ex(3) = -(                                              4*prestencil_ex(2)+ 8*prestencil_ex(3))/256
    stencilz_ex(4) =  (                                                                  1*prestencil_ex(3))/256

    stencilz_by(0) =  (256+128*prestencil_by(0)+96*prestencil_by(1)+80*prestencil_by(2)+70*prestencil_by(3))/256
    stencilz_by(1) = -(     64*prestencil_by(0)+64*prestencil_by(1)+60*prestencil_by(2)+56*prestencil_by(3))/256
    stencilz_by(2) =  (                         16*prestencil_by(1)+24*prestencil_by(2)+28*prestencil_by(3))/256
    stencilz_by(3) = -(                                              4*prestencil_by(2)+ 8*prestencil_by(3))/256
    stencilz_by(4) =  (                                                                  1*prestencil_by(3))/256
  ELSE
    WRITE(*,*) "ERROR: NCI corrector requires l_lower_order_in_v=1, i.e., Galerkin scheme"
    STOP
  ENDIF
  
END SUBROUTINE init_godfrey_filter_coeffs

! ________________________________________________________________________________________
!> @brief
!> Apply a stencil along z for a 2d array
!> 
!> @param[in] lo, hi     1d, 2-element arrays: lowest and highest indices between which the 
!>                       filter is applied
!> @param[inout] fou     2d (x-z) array: Output filtered field
!> @param[in] olo, ohi   1d, 2-element arrays: lower and upper bounds for the filtered 
!>                       output field
!> @param[in] fin        2d (x-z) array: Input field to filter
!> @param[in] ilo, ihi   1d, 2-element arrays: lower and upper bounds for the input field 
!>                       to filter
!> @param[in] sten       1d stencil along z
!> @param[in] nsten      stencil order
!>
!> @date
!> Creation 2018
! ________________________________________________________________________________________
SUBROUTINE apply_godfrey_filter_z_2d (lo, hi, fou, olo, ohi, fin, ilo, ihi, sten, nsten) bind(c,name='apply_godfrey_filter_z_2d')
  INTEGER, INTENT(in) :: lo(2), hi(2), olo(2), ohi(2), ilo(2), ihi(2), nsten
  REAL(NUM), INTENT(INOUT) :: fou(olo(1):ohi(1),olo(2):ohi(2))
  REAL(NUM), INTENT(IN   ) :: fin(ilo(1):ihi(1),ilo(2):ihi(2))
  REAL(NUM), INTENT(IN   ) :: sten(0:nsten-1)
  INTEGER :: i,k,ks

  DO    k = lo(2), hi(2)
     DO i = lo(1), hi(1)
        fou(i,k) = sten(0)*fin(i,k)
     END DO
  END DO

  DO ks = 1, nsten-1
     DO    k = lo(2), hi(2)
        DO i = lo(1), hi(1)
           fou(i,k) = fou(i,k) + sten(ks)*(fin(i,k-ks)+fin(i,k+ks))
        END DO
     END DO
  END DO
END SUBROUTINE apply_godfrey_filter_z_2d

! ________________________________________________________________________________________
!> @brief
!> Apply a stencil along z for a 3d array
!> 
!> @param[in] lo, hi     1d, 3-element arrays: lowest and highest indices between which the 
!>                       filter is applied
!> @param[inout] fou     3d (x-y-z) array: Output filtered field
!> @param[in] olo, ohi   1d, 3-element arrays: lower and upper bounds for the filtered 
!>                       output field
!> @param[in] fin        3d (x-y-z) array: Input field to filter
!> @param[in] ilo, ihi   1d, 3-element arrays: lower and upper bounds for the input field 
!>                       to filter
!> @param[in] sten       1d stencil along z
!> @param[in] nsten      stencil order
!>
!> @date
!> Creation 2018
! ________________________________________________________________________________________
SUBROUTINE apply_godfrey_filter_z_3d (lo, hi, fou, olo, ohi, fin, ilo, ihi, sten, nsten) bind(c,name='apply_godfrey_filter_z_3d')
  INTEGER, INTENT(in) :: lo(3), hi(3), olo(3), ohi(3), ilo(3), ihi(3), nsten
  REAL(NUM), INTENT(INOUT) :: fou(olo(1):ohi(1),olo(2):ohi(2),olo(3):ohi(3))
  REAL(NUM), INTENT(IN   ) :: fin(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3))
  REAL(NUM), INTENT(IN   ) :: sten(0:nsten-1)
  INTEGER :: i,j,k,ks

  DO       k = lo(3), hi(3)
     DO    j = lo(2), hi(2)
        DO i = lo(1), hi(1)
           fou(i,j,k) = sten(0)*fin(i,j,k)
        END DO
     END DO
  END DO
  
  DO ks = 1, nsten-1
     DO       k = lo(3), hi(3)
        DO    j = lo(2), hi(2)
           DO i = lo(1), hi(1)
              fou(i,j,k) = fou(i,j,k) + sten(ks)*(fin(i,j,k-ks)+fin(i,j,k+ks))
           END DO
        END DO
     END DO
  END DO
END SUBROUTINE apply_godfrey_filter_z_3d

END MODULE godfrey_filter_coeffs
