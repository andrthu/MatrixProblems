set(CMAKE_C_COMPILER "/usr/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "11.4.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "17")
set(CMAKE_C_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert;c_std_17;c_std_23")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")
set(CMAKE_C17_COMPILE_FEATURES "c_std_17")
set(CMAKE_C23_COMPILE_FEATURES "c_std_23")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_C_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "/usr/bin/gcc-ar-11")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/usr/bin/gcc-ranlib-11")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)

set(CMAKE_C_COMPILER_ENV_VAR "CC")

set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_C_LIBRARY_ARCHITECTURE "aarch64-linux-gnu")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "aarch64-linux-gnu")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/cm/shared/apps/slurm/current/include;/home/andreast/modules/arm/apps/flow/202309a/rpp-amgcpr-ghost-last/include;/home/andreast/modules/arm/apps/zoltan/202309a/12.16/include;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/suitesparse-32-5.12.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/mpfr-4.0.2/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/gmp-6.3.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/parmetis-32-4.0.3/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/metis-32-5.1.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/openblas-0.3.21/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/boost-1.73.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/openmpi-4.1.4/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/util-linux-2.34/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libfabric-1.17.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/ucx-1.12.1/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/rdma-core-44.0/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libnl-3.2.25/include/libnl3;/cm/shared/ex3-modules/202309a/huaq/pkgs/numactl-2.0.13/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libevent-2.1.12-stable/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/knem-1.1.4/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/hwloc-2.7.1/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/libxml2-2.9.12/include;/cm/shared/ex3-modules/202309a/huaq/pkgs/xz-5.2.7/include;/usr/lib/gcc/aarch64-linux-gnu/11/include;/usr/local/include;/usr/include/aarch64-linux-gnu;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/home/andreast/modules/arm/apps/flow/202309a/rpp-amgcpr-ghost-last/lib;/home/andreast/modules/arm/apps/zoltan/202309a/12.16/lib;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/suitesparse-32-5.12.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/mpfr-4.0.2/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/gmp-6.3.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/parmetis-32-4.0.3/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/metis-32-5.1.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/openblas-0.3.21/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/boost-1.73.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/openmpi-4.1.4/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/util-linux-2.34/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libfabric-1.17.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/ucx-1.12.1/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/rdma-core-44.0/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libnl-3.2.25/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/numactl-2.0.13/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libevent-2.1.12-stable/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/knem-1.1.4/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/hwloc-2.7.1/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/libxml2-2.9.12/lib;/cm/shared/ex3-modules/202309a/huaq/pkgs/xz-5.2.7/lib;/usr/lib/gcc/aarch64-linux-gnu/11;/usr/lib/aarch64-linux-gnu;/usr/lib;/lib/aarch64-linux-gnu;/lib;/cm/shared/apps/slurm/current/lib64/slurm;/cm/shared/apps/slurm/current/lib64")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
