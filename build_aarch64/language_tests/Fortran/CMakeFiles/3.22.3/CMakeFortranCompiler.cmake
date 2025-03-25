set(CMAKE_Fortran_COMPILER "/usr/bin/f95")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "11.4.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/usr/bin/gcc-ar-11")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/usr/bin/gcc-ranlib-11")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

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
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "aarch64-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "aarch64-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/usr/lib/gcc/aarch64-linux-gnu/11/finclude;/cm/shared/apps/slurm/current/include;/home/andreast/modules/arm/apps/flow/202309a/rpp-amgcpr-ghost-last/include;/home/andreast/modules/arm/apps/zoltan/202309a/12.16/include;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/suitesparse-32-5.12.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/mpfr-4.0.2/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/gmp-6.3.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/parmetis-32-4.0.3/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/metis-32-5.1.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/openblas-0.3.21/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/boost-1.73.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/openmpi-4.1.4/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/util-linux-2.34/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libfabric-1.17.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/ucx-1.12.1/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/rdma-core-44.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libnl-3.2.25/include/libnl3;/cm/shared/ex3-modules/202309a/huaq/pkgs/numactl-2.0.13/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libevent-2.1.12-stable/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/knem-1.1.4/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/hwloc-2.7.1/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libxml2-2.9.12/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/xz-5.2.7/include;/usr/lib/gcc/aarch64-linux-gnu/11/include;/usr/local/include;/usr/include/aarch64-linux-gnu;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/home/andreast/modules/arm/apps/flow/202309a/rpp-amgcpr-ghost-last/lib;/home/andreast/modules/arm/apps/zoltan/202309a/12.16/lib;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/suitesparse-32-5.12.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/mpfr-4.0.2/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/gmp-6.3.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/parmetis-32-4.0.3/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/metis-32-5.1.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/openblas-0.3.21/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/boost-1.73.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/openmpi-4.1.4/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/util-linux-2.34/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libfabric-1.17.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/ucx-1.12.1/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/rdma-core-44.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libnl-3.2.25/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/numactl-2.0.13/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libevent-2.1.12-stable/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/knem-1.1.4/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/hwloc-2.7.1/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libxml2-2.9.12/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/xz-5.2.7/lib;/usr/lib/gcc/aarch64-linux-gnu/11;/usr/lib/aarch64-linux-gnu;/usr/lib;/lib/aarch64-linux-gnu;/lib;/cm/shared/apps/slurm/current/lib64/slurm;/cm/shared/apps/slurm/current/lib64")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
