# Install script for directory: /home/erikhs/MatrixProblems

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  set(CMAKE_MODULE_PATH /home/erikhs/MatrixProblems/cmake/modules;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/share/dune/cmake/modules;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/share/dune/cmake/modules;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/share/dune/cmake/modules;/home/andreast/modules/arm/apps/dune/202309a/mod-2.7.1/share/dune/cmake/modules)
              set(DUNE_PYTHON_WHEELHOUSE /usr/local/share/dune/wheelhouse)
              include(DuneExecuteProcess)
              dune_execute_process(COMMAND "/cm/shared/ex3-modules/202309a/huaq/pkgs/cmake-3.22.3/bin/cmake" --build . --target install_python --config $<CONFIG>)
              
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/dunecontrol/MatrixProblems" TYPE FILE FILES "/home/erikhs/MatrixProblems/dune.module")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/MatrixProblems" TYPE FILE FILES
    "/home/erikhs/MatrixProblems/build_aarch64/cmake/pkg/MatrixProblems-config.cmake"
    "/home/erikhs/MatrixProblems/build_aarch64/MatrixProblems-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/MatrixProblems" TYPE FILE FILES "/home/erikhs/MatrixProblems/config.h.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/erikhs/MatrixProblems/build_aarch64/MatrixProblems.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/erikhs/MatrixProblems/build_aarch64/src/cmake_install.cmake")
  include("/home/erikhs/MatrixProblems/build_aarch64/m4/cmake_install.cmake")
  include("/home/erikhs/MatrixProblems/build_aarch64/dune/cmake_install.cmake")
  include("/home/erikhs/MatrixProblems/build_aarch64/doc/cmake_install.cmake")
  include("/home/erikhs/MatrixProblems/build_aarch64/cmake/modules/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/erikhs/MatrixProblems/build_aarch64/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
