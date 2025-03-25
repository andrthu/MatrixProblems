set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "11.4.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/bin/gcc-ar-11")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/bin/gcc-ranlib-11")
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
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/home/andreast/modules/flow/202309a/rpp-dune-amgcpr-branch/lib;/home/andreast/modules/zoltan/202309a/12.16/lib;/home/andreast/modules/dune/202309a/mod-2.7.1/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/suitesparse-32-5.12.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/mpfr-4.0.2/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/gmp-6.1.2/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/parmetis-32-4.0.3/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/metis-32-5.1.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/boost-1.73.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/xz-5.2.7/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/openblas-0.3.21/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/openmpi-4.1.4/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/util-linux-2.34/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/libfabric-1.17.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/ucx-1.12.1/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/rdma-core-44.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/libnl-3.2.25/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/numactl-2.0.13/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/libevent-2.1.12-stable/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/knem-1.1.4/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/hwloc-2.7.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/fftw-3.3.8/lib;/cm/shared/ex3-modules/0.6.1/pkgs/lapack-3.9.0/lib;/cm/shared/ex3-modules/0.6.1/pkgs/openblas-0.3.12/lib;/cm/shared/ex3-modules/0.6.1/pkgs/gcc-8.4.0/lib;/cm/shared/ex3-modules/0.6.1/pkgs/freetype-2.10.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/libpng-1.6.37/lib;/cm/shared/ex3-modules/0.6.1/pkgs/python-3.7.4/lib;/cm/shared/ex3-modules/0.6.1/pkgs/sqlite-3.31.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/libffi-3.2.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/readline-8.0/lib;/cm/shared/ex3-modules/0.6.1/pkgs/ncurses-6.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/xz-5.2.5/lib;/cm/shared/ex3-modules/0.6.1/pkgs/bzip2-1.0.8/lib;/cm/shared/ex3-modules/0.6.1/pkgs/openssl-1.1.1c/lib;/cm/shared/apps/slurm/current/lib;/usr/lib/gcc/x86_64-linux-gnu/11;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib;/cm/shared/ex3-modules/0.6.1/pkgs/gcc-8.4.0/lib64;/cm/shared/apps/slurm/current/lib64/slurm;/cm/shared/apps/slurm/current/lib64")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
