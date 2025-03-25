set(CMAKE_C_COMPILER "/usr/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "11.4.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

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
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/cm/shared/apps/slurm/current/include;/home/andreast/modules/flow/202309a/rpp-dune-amgcpr-branch/include;/home/andreast/modules/zoltan/202309a/12.16/include;/home/andreast/modules/dune/202309a/mod-2.7.1/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/suitesparse-32-5.12.0/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/mpfr-4.0.2/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/gmp-6.1.2/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/parmetis-32-4.0.3/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/metis-32-5.1.0/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/boost-1.73.0/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/xz-5.2.7/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/openblas-0.3.21/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/openmpi-4.1.4/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/util-linux-2.34/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/libfabric-1.17.0/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/ucx-1.12.1/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/rdma-core-44.0/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/libnl-3.2.25/include/libnl3;/cm/shared/ex3-modules/202309a/milanq/pkgs/numactl-2.0.13/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/libevent-2.1.12-stable/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/knem-1.1.4/include;/cm/shared/ex3-modules/202309a/milanq/pkgs/hwloc-2.7.1/include;/cm/shared/ex3-modules/0.6.1/pkgs/fftw-3.3.8/include;/cm/shared/ex3-modules/0.6.1/pkgs/lapack-3.9.0/include;/cm/shared/ex3-modules/0.6.1/pkgs/openblas-0.3.12/include;/cm/shared/ex3-modules/0.6.1/pkgs/freetype-2.10.1/include;/cm/shared/ex3-modules/0.6.1/pkgs/libpng-1.6.37/include;/cm/shared/ex3-modules/0.6.1/pkgs/python-3.7.4/include/python3.7m;/cm/shared/ex3-modules/0.6.1/pkgs/sqlite-3.31.1/include;/cm/shared/ex3-modules/0.6.1/pkgs/libffi-3.2.1/include;/cm/shared/ex3-modules/0.6.1/pkgs/readline-8.0/include;/cm/shared/ex3-modules/0.6.1/pkgs/ncurses-6.1/include;/cm/shared/ex3-modules/0.6.1/pkgs/xz-5.2.5/include;/cm/shared/ex3-modules/0.6.1/pkgs/bzip2-1.0.8/include;/cm/shared/ex3-modules/0.6.1/pkgs/openssl-1.1.1c/include;/usr/lib/gcc/x86_64-linux-gnu/11/include;/usr/local/include;/usr/include/x86_64-linux-gnu;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/home/andreast/modules/flow/202309a/rpp-dune-amgcpr-branch/lib;/home/andreast/modules/zoltan/202309a/12.16/lib;/home/andreast/modules/dune/202309a/mod-2.7.1/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/suitesparse-32-5.12.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/mpfr-4.0.2/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/gmp-6.1.2/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/parmetis-32-4.0.3/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/metis-32-5.1.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/boost-1.73.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/xz-5.2.7/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/openblas-0.3.21/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/openmpi-4.1.4/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/util-linux-2.34/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/libfabric-1.17.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/ucx-1.12.1/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/rdma-core-44.0/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/libnl-3.2.25/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/numactl-2.0.13/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/libevent-2.1.12-stable/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/knem-1.1.4/lib;/cm/shared/ex3-modules/202309a/milanq/pkgs/hwloc-2.7.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/fftw-3.3.8/lib;/cm/shared/ex3-modules/0.6.1/pkgs/lapack-3.9.0/lib;/cm/shared/ex3-modules/0.6.1/pkgs/openblas-0.3.12/lib;/cm/shared/ex3-modules/0.6.1/pkgs/gcc-8.4.0/lib;/cm/shared/ex3-modules/0.6.1/pkgs/freetype-2.10.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/libpng-1.6.37/lib;/cm/shared/ex3-modules/0.6.1/pkgs/python-3.7.4/lib;/cm/shared/ex3-modules/0.6.1/pkgs/sqlite-3.31.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/libffi-3.2.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/readline-8.0/lib;/cm/shared/ex3-modules/0.6.1/pkgs/ncurses-6.1/lib;/cm/shared/ex3-modules/0.6.1/pkgs/xz-5.2.5/lib;/cm/shared/ex3-modules/0.6.1/pkgs/bzip2-1.0.8/lib;/cm/shared/ex3-modules/0.6.1/pkgs/openssl-1.1.1c/lib;/cm/shared/apps/slurm/current/lib;/usr/lib/gcc/x86_64-linux-gnu/11;/usr/lib/x86_64-linux-gnu;/usr/lib;/lib/x86_64-linux-gnu;/lib;/cm/shared/ex3-modules/0.6.1/pkgs/gcc-8.4.0/lib64;/cm/shared/apps/slurm/current/lib64/slurm;/cm/shared/apps/slurm/current/lib64")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
