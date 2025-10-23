set(CMAKE_MODULE_PATH /home/andreast/modules/flow/202309a/nov2024Alf/share/opm/cmake/Modules)
#cmake_policy(SET CMP0028 OLD)

cmake_policy(SET CMP0060 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0028 OLD)

find_package(Zoltan CONFIG REQUIRED)
find_package(opm-common CONFIG REQUIRED)
find_package(opm-grid CONFIG REQUIRED)
find_package(ParMETIS REQUIRED)


find_package(opm-simulators CONFIG REQUIRED)

add_executable("MatrixProblems" MatrixProblems.cc)

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

add_executable("wellJsonSolve" 
  wellJsonSolve.cpp
  io/inputRead.hpp
  io/readWells.hpp
  io/standardWell.hpp
  io/multiSegmentWell.hpp
  io/dictRead.hpp
  solveParallel/ghostLastOperations.hpp
  solveParallel/operatorWithWell.hpp
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
  generalDim/wellJsonSolveGeneral.hpp)

add_executable("eclipseIO"
  eclipseIO.cpp
  )

set(opm_dep
  jsonSolve
  jsonSolveCO2
  multJsonSolve
  multJsonSolveCO2
  singleHierarchyMultJsonSolve
  singleHierarchyMultJsonSolveCO2
  wellJsonSolve
  eclipseIO
  )


INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

foreach(exe ${opm_dep})

  target_link_libraries(${exe}  ${opm-simulators_LIBRARIES} )
  target_link_libraries(${exe}  ${opm-grid_LIBRARIES} )

  target_link_libraries(${exe}  ${opm-common_LIBRARIES} )

  target_link_libraries(${exe} ${BOOST_LIBRARYDIR}/libboost_program_options.so )
  target_link_libraries(${exe} ${BOOST_LIBRARYDIR}/libboost_filesystem.so )

  
endforeach()
