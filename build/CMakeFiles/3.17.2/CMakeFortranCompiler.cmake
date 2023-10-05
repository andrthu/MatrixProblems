set(CMAKE_Fortran_COMPILER "/cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/gfortran-10.1")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "10.1.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/gcc-ar-10.1")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/gcc-ranlib-10.1")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/lib64;/cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/lib/gcc/x86_64-pc-linux-gnu/10.1.0;/lib/x86_64-linux-gnu;/lib64;/usr/lib/x86_64-linux-gnu;/usr/lib64;/home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib;/cm/shared/ex3-modules/0.5.0/pkgs/scotch-6.0.7/lib;/cm/shared/ex3-modules/0.4.0/apps/python/3.7.4/lib;/cm/shared/ex3-modules/0.4.0/apps/sqlite/3.31.1/lib;/cm/shared/ex3-modules/0.4.0/apps/readline/8.0/lib;/cm/shared/ex3-modules/0.4.0/apps/ncurses/6.1/lib;/cm/shared/ex3-modules/0.4.0/apps/xz/5.2.5/lib;/cm/shared/ex3-modules/0.4.0/apps/bzip2/1.0.8/lib;/cm/shared/ex3-modules/0.5.0/pkgs/fftw-3.3.8/lib;/cm/shared/ex3-modules/0.5.0/pkgs/lapack-3.9.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/cblas/lib;/cm/shared/ex3-modules/0.5.0/pkgs/netlib-blas-3.8.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/freetype-2.10.1/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libpng-1.6.37/lib;/cm/shared/ex3-modules/0.5.0/pkgs/python-3.7.4/lib;/cm/shared/ex3-modules/0.5.0/pkgs/sqlite-3.31.1/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libffi-3.2.1/lib;/cm/shared/ex3-modules/0.5.0/pkgs/bzip2-1.0.8/lib;/cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/xz-5.2.5/lib;/home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib;/home/andreast/modules/parmetis/0.5.0/4.0.3/lib;/cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib;/cm/shared/ex3-modules/0.5.0/pkgs/metis-5.1.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib;/home/andreast/modules/zoltan/0.5.0/12.16/lib;/cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libfabric-1.11.1/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libevent-2.1.11-stable/lib;/cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/mpc-1.1.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/mpfr-4.0.2/lib;/cm/shared/ex3-modules/0.5.0/pkgs/gmp-6.1.2/lib;/cm/shared/ex3-modules/0.5.0/pkgs/slurm-20.02.5/lib;/cm/shared/ex3-modules/0.5.0/pkgs/readline-8.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/ncurses-6.1/lib;/cm/shared/ex3-modules/0.5.0/pkgs/curl-7.69.1/lib;/cm/shared/ex3-modules/0.5.0/pkgs/munge-0.5.13/lib;/cm/shared/ex3-modules/0.5.0/pkgs/openssl-1.1.1c/lib;/cm/shared/ex3-modules/0.5.0/pkgs/freeipmi-1.6.6/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libgcrypt-1.8.7/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libgpg-error-1.39/lib;/cm/shared/ex3-modules/0.5.0/pkgs/ucx-1.9.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/rdma-core-31.0/lib;/cm/shared/ex3-modules/0.5.0/pkgs/libnl-3.2.25/lib;/cm/shared/ex3-modules/0.5.0/pkgs/numactl-2.0.13/lib;/cm/shared/ex3-modules/0.5.0/pkgs/knem-1.1.4/lib;/cm/shared/ex3-modules/0.5.0/pkgs/hwloc-2.0.4/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
