set(CMAKE_MODULE_PATH /home/andreast/modules/flow/202309a/ss/master2024.oct/share/opm/cmake/Modules)
#cmake_policy(SET CMP0028 OLD)

cmake_policy(SET CMP0060 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0028 OLD)

find_package(Zoltan CONFIG REQUIRED)
find_package(opm-common CONFIG REQUIRED)
find_package(opm-grid CONFIG REQUIRED)
#find_package(METIS REQUIRED)
find_package(ParMETIS REQUIRED)
#add_library(METIS::METIS UNKNOWN IMPORTED)
#set_target_properties(METIS::METIS PROPERTIES
#  IMPORTED_LOCATION /cm/shared/ex3-modules/202309a/milanq/pkgs/metis-32-5.1.0/lib/libmetis.so
#  INTERFACE_INCLUDE_DIRECTORIES /cm/shared/ex3-modules/202309a/milanq/pkgs/metis-32-5.1.0/include
#  )##
#add_library(ParMETIS::ParMETIS UNKNOWN IMPORTED)
#set_target_properties(ParMETIS::ParMETIS PROPERTIES
#    IMPORTED_LOCATION /cm/shared/ex3-modules/202309a/milanq/pkgs/parmetis-32-4.0.3/lib/libparmetis.so
#    INTERFACE_INCLUDE_DIRECTORIES /cm/shared/ex3-modules/202309a/milanq/pkgs/parmetis-32-4.0.3/include
#    )

find_package(OpenMP REQUIRED)


find_package(opm-simulators CONFIG REQUIRED)

add_executable("MatrixProblems" MatrixProblems.cc)

#add_executable("twoLevelMethodTest" 
#  twoLevelMethodTest.cpp
#  io/inputRead.hpp
#  io/dictRead.hpp
#  solveParallel/ghostLastOperations.hpp
#  partition/graphFunctions.hpp
#  partition/overlapCreation.hpp
#  solveParallel/reorder.hpp
#  solveParallel/buildLocalMatrix.hpp
#  solveParallel/nonQuadraticLocalMatrix.hpp
#  partition/setIndexSet.hpp
#  partition/evaluatePartition.hpp
#  partition/meassureCommunication.hpp
#  partition/meassureOperationsMinLoop.hpp
#  partition/flexibleSolverTimer.hpp
#  generalDim/twoLevelMethod.hpp)

add_executable("jsonSolve" 
  jsonSolve.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp)

add_executable("jsonSolveCO2" 
  jsonSolveCO2.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp)

add_executable("hirTest" 
  hirTest.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp
  generalDim/hirTest.hpp)

add_executable("hirTest1dim" 
  hirTest1dim.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp
  generalDim/hirTest.hpp)

add_executable("listTransInfo"
  listTransInfo.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp
  generalDim/listTransInfo.hpp)

add_executable("listTransInfoCo2"
  listTransInfoCo2.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp
  generalDim/listTransInfo.hpp)

add_executable("zoltanInfo"
  zoltanInfo.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp
  generalDim/zoltanInfo.hpp)

add_executable("zoltanInfoCo2"
  zoltanInfoCo2.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/jsonSolveGeneral.hpp
  generalDim/zoltanInfo.hpp)

add_executable("multJsonSolve" 
  multJsonSolve.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/multJsonSolveGeneral.hpp)

add_executable("multJsonSolveCO2" 
  multJsonSolveCO2.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/multJsonSolveGeneral.hpp)

add_executable("singleHierarchyMultJsonSolve" 
  singleHierarchyMultJsonSolve.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/singleHierarchymultJsonSolveGeneral.hpp)

add_executable("singleHierarchyMultJsonSolveCO2" 
  singleHierarchyMultJsonSolveCO2.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/singleHierarchymultJsonSolveGeneral.hpp)

add_executable("transmissibilityHirPart" 
  transmissibilityHirPart.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/transmissibilityHirPart.hpp)

add_executable("transmissibilityHirPartCo2"
  transmissibilityHirPartCo2.cpp
  io/inputRead.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  partition/graphFunctions.hpp
  partition/overlapCreation.hpp
  solveParallel/reorder.hpp
  solveParallel/buildLocalMatrix.hpp
  solveParallel/nonQuadraticLocalMatrix.hpp
  partition/setIndexSet.hpp
  partition/evaluatePartition.hpp
  partition/meassureCommunication.hpp
  partition/meassureOperationsMinLoop.hpp
  partition/flexibleSolverTimer.hpp
  generalDim/transmissibilityHirPart.hpp)

add_executable("eclipseIO"
  eclipseIO.cpp
  )

set(opm_dep
  #twoLevelMethodTest
  jsonSolve
  jsonSolveCO2
  hirTest
  hirTest1dim
  multJsonSolve
  multJsonSolveCO2
  singleHierarchyMultJsonSolve
  singleHierarchyMultJsonSolveCO2
  eclipseIO
  listTransInfo
  listTransInfoCo2
  zoltanInfo
  zoltanInfoCo2
  transmissibilityHirPart
  transmissibilityHirPartCo2
  )


INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

foreach(exe ${opm_dep})

  target_link_libraries(${exe}  ${opm-simulators_LIBRARIES} )
  target_link_libraries(${exe}  ${opm-grid_LIBRARIES} )
  target_link_libraries(${exe}  ${opm-common_LIBRARIES} )

  target_link_libraries(${exe} ${BOOST_LIBRARYDIR}/libboost_program_options.so )
  target_link_libraries(${exe} ${BOOST_LIBRARYDIR}/libboost_filesystem.so )
  #target_link_libraries(${exe} PRIVATE OpenMP::OpenMP_CXX)
  
endforeach()
