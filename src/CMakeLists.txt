add_executable("MatrixProblems" MatrixProblems.cc)
cmake_policy(SET CMP0060 NEW)
#cmake_policy(SET CMP0028 OLD)
#cmake_policy(SET CMP0074 OLD)

find_package(Zoltan CONFIG REQUIRED)
find_package(opm-common CONFIG REQUIRED)
find_package(opm-grid CONFIG REQUIRED)
find_package(opm-material CONFIG REQUIRED)
#find_package(ewoms CONFIG REQUIRED)
find_package(opm-simulators CONFIG REQUIRED)
#find_package(METIS REQUIRED)
#find_package(ParMETIS REQUIRED)

#add_library(METIS::METIS UNKNOWN IMPORTED)
#set_target_properties(METIS::METIS PROPERTIES
#  IMPORTED_LOCATION ${METIS_LIBRARY}
#  INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIR}
#  INTERFACE_COMPILE_DEFINITIONS METIS_API_VERSION=${METIS_API_VERSION}
#  )

#dune_register_package_flags(LIBRARIES METIS::METIS)

#target_link_dune_default_libraries("MatrixProblems")

add_executable("readMatrixProduceInfo" 
   readMatrixProduceInfo.cpp
   io/inputRead.hpp
   io/dictRead.hpp
   solveParallel/ghostLastOperations.hpp
   partition/graphFunctions.hpp
   partition/overlapCreation.hpp
   solveParallel/reorder.hpp
   solveParallel/buildLocalMatrix.hpp
   partition/setIndexSet.hpp
   partition/evaluatePartition.hpp
   solveParallel/solveSystem.hpp
   solveParallel/compareSolutions.hpp
   partition/meassureCommunication.hpp
   partition/timeMultipleTimes.hpp)

target_link_dune_default_libraries("readMatrixProduceInfo")
if(TRUE)
add_executable("wellReadAndSolve" 
   wellReadAndSolve.cpp
   io/inputRead.hpp
   io/wells.hpp
   io/dictRead.hpp
   solveParallel/ghostLastOperations.hpp
   partition/graphFunctions.hpp
   partition/overlapCreation.hpp
   solveParallel/reorder.hpp
   solveParallel/buildLocalMatrix.hpp
   partition/setIndexSet.hpp
   partition/evaluatePartition.hpp
   solveParallel/solveSystem.hpp
   solveParallel/wellSolveSystem.hpp
   solveParallel/compareSolutions.hpp
   partition/meassureCommunication.hpp
   partition/timeMultipleTimes.hpp)

target_link_dune_default_libraries("wellReadAndSolve")

add_executable("meassureComponents" 
   meassureComponents.cpp
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
   partition/meassureOperationsMinLoop.hpp)

add_executable("amgMes" 
   amg_meassure.cpp
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
   partition/meassureOperationsMinLoop.hpp)

add_executable("flexibleSolver" 
   flexibleSolver.cpp
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
   partition/flexibleSolverTimer.hpp)


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
  generalDim/GhostLastTest.hpp
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
  generalDim/GhostLastTest.hpp
  generalDim/jsonSolveGeneral.hpp)

add_executable("oneIterTest" 
   oneIterTest.cpp
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
   partition/flexibleSolverTimer.hpp)

add_executable("amgCprTest"
  amgCprTest.cpp
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
  solveParallel/amgSetup.hpp
  generalDim/amgCpr.hpp)

add_executable("amgCprCo2Test"
  amgCprCo2Test.cpp
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
  solveParallel/amgSetup.hpp)

add_executable("amgCprGLTest"
  amgCprGLTest.cpp
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
  solveParallel/amgSetup.hpp
  solveParallel/copyILU0.hpp
  generalDim/GhostLastTest.hpp)


add_executable("updateILU"
  updateILU.cpp
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
  solveParallel/amgSetup.hpp
  solveParallel/copyILU0.hpp
  generalDim/updateTest.hpp)


add_executable("amgCprGLCo2Test"
  amgCprGLCo2Test.cpp
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
  solveParallel/amgSetup.hpp
  solveParallel/copyILU0.hpp
  generalDim/GhostLastTest.hpp)

if(TRUE)
add_executable("twoLevelMethodTest"
  twoLevelMethodTest.cpp
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
  solveParallel/amgSetup.hpp
  generalDim/twoLevelMethod.hpp)
#endif()

add_executable("tunedAmg"
  tunedAmg.cpp
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
  solveParallel/amgSetup.hpp
  generalDim/onlyAMG.hpp)

add_executable("tunedAmgCo2"
  tunedAmgCo2.cpp
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
  solveParallel/amgSetup.hpp
  generalDim/onlyAMG.hpp)

#if(FALSE)
add_executable("twoLevelMethodTestCO2"
  twoLevelMethodTestCO2.cpp
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
  solveParallel/amgSetup.hpp
  generalDim/twoLevelMethod.hpp)
endif()

add_executable("co2FS" 
   co2FS.cpp
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
   partition/flexibleSolverTimer.hpp)

add_executable("doTwoPressureIndex" 
   doTwoPressureIndex.cpp
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
   partition/flexibleSolverTimer.hpp)

add_executable("betterIOTest" 
   betterIOTest.cpp
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
   partition/flexibleSolverTimer.hpp)

add_executable("betterIOTest2" 
  betterIOTest2.cpp
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
  solveParallel/copyILU0.hpp)

add_executable("blockMatComp" 
  blockMatComp.cpp
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
  partition/flexibleSolverTimer.hpp)

 
add_executable("eclipseIO"
  eclipseIO.cpp
  )

target_link_dune_default_libraries("meassureComponents")

endif()
add_executable("justPartition" 
   justPartition.cpp
   io/inputRead.hpp
   io/dictRead.hpp
   partition/graphFunctions.hpp
   partition/overlapCreation.hpp
   partition/setIndexSet.hpp
   partition/evaluatePartition.hpp)
   


set(my_EXECUTABLES
  readMatrixProduceInfo
  justPartition
  )
if (TRUE)
set(opm_dep
  wellReadAndSolve
  meassureComponents
  amgMes
  flexibleSolver
  jsonSolve
  jsonSolveCO2
  oneIterTest
  amgCprTest
  amgCprCo2Test
  amgCprGLTest
  amgCprGLCo2Test
  twoLevelMethodTest
  twoLevelMethodTestCO2
  tunedAmg
  tunedAmgCo2
  co2FS
  updateILU
  doTwoPressureIndex
  betterIOTest
  betterIOTest2
  blockMatComp
  eclipseIO
  )

endif()

foreach(exe ${my_EXECUTABLES})

  target_link_libraries(${exe}  ${Zoltan_LIBRARIES} )

  target_link_libraries(${exe}  ${MPI_LIBRARIES} )
  target_link_libraries(${exe}  /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so )
endforeach()


foreach(exe ${opm_dep})

  target_link_libraries(${exe}  ${opm-simulators_LIBRARIES} )
  target_link_libraries(${exe}  ${opm-grid_LIBRARIES} )
  target_link_libraries(${exe}  ${opm-material_LIBRARIES} )
  target_link_libraries(${exe}  ${opm-common_LIBRARIES} )

  target_link_libraries(${exe}  /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so )
  
endforeach()
