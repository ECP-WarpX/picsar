# **Installation of the PICSAR Python full mode**

In order to use highly scalable hybrid PSATD with WARP through PICSAR library, 
you will need to install FFTW_MPI and P3DFFT library.


I-Install FFTW_MPI Library:
To install FFTW_MPI you can follow the instructions here :
http://www.fftw.org/fftw2_doc/fftw_6.html
Use this configuration (for gnu compilation):

./configure --enable-mpi --enable-openmp --enable-shared --enable-gnu --prefix=/install_dir_fftw
Then
-make
-make install
and
-export LD_LIBRARY_FILE=LD_LIBRARY_FILE:install_dir_fftw/lib
II-P3DFFT library can be downloaded here:



"https://www.p3dfft.net/"




To get good performances, you need to install p3dfft with this configuration:
-./configure --enable-stride1 --enable-openmp --enable-measure --enable-fftw --enable-gnu --with-fftw=$FFTW_DIR --prefix= path_to_install_directory

then:
-make
-make install

To use p3dfft along with warp though picsar, p3dfft needs to be linked dynamically to picsar
But the installation instructions above will generate a static library.
What you can do is to install it manually.

The files needed for compilation are : 
   fft_init.F90, fft_exec.F90, wrap.F90, fft_spec.F90, module.F90
Check the installation logs (when you run make) and look for the compilation lines for these files.
They should look like this:
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include   -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o fft_spec.o fft_spec.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include   -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o module.o module.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include   -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o fft_init.o fft_init.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include   -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o fft_exec.o fft_exec.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include   -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o wrap.o wrap.F90

enter build/ directory and recompile these files adding -fPIC flag

-cd build/
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include -fPIC  -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o fft_spec.o fft_spec.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include -fPIC  -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o module.o module.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include -fPIC  -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o fft_init.o fft_init.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include -fPIC  -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o fft_exec.o fft_exec.F90
-mpif90 -DHAVE_CONFIG_H -I. -I.. -I/home/install_fftw/include -fPIC  -DGNU -DOPENMP -fopenmp -DMEASURE -DSTRIDE1 -DFFTW -g -O2 -c -o wrap.o wrap.F90

and then generate a dynamic link to p3dfft

-mpif90 -shared  fft_spec.o module.o fft_init.o fft_exec.o wrap.o -o libp3dfft.so
this should generate a libp3dfft.so
Copy libp3dfft.so to the path_to_install_directory/lib
and 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_install_directory/lib (you may want to add this line into your .bashrc)





## Editing configure file:

In configure file line 36: 
Set 
-fftw_lib="/path_to_fftw_directory/lib"
-p3dfft_lib="/path_to_p3dfft_directory/lib"
-fftw_link ="-I /path_to_fftw_directory/include"
-p3dfft_link = "-I path_to_p3dfft_directory/include" 
Finally type:

python configure --full_pxr True --other_flags(compiler - architecture) ...


