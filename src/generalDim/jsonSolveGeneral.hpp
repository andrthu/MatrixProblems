/*
  Copyright 2023 Andreas Thune.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_JSONSOLVEGENERAL_HEADER_INCLUDED
#define OPM_JSONSOLVEGENERAL_HEADER_INCLUDED

#endif // OPM_JSONSOLVEGENERAL_HEADER_INCLUDED

template<class Mat, class Vec>
void gen_dim_jsonSolve(int argc, char** argv)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;                 // solveParallel/ghostLastOperations.hpp
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    const auto block_size = Vec::block_type::dimension;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    //findZeroDiag(A_loc, rhs_loc);
    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);

    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
    	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);//0);//10);
    }
    prm_json.put("preconditioner.verbosity", 10);//0);//10);
    prm_json.put("verbosity", 2);//0);//2);
    if (pc_Type == "cpr") {
        prm_json.put("preconditioner.coarsesolver.preconditioner.maxlevel", 15);
        prm_json.put("preconditioner.coarsesolver.maxiter", 1);
    } else if (pc_Type == "amg") {
        prm_json.put("preconditioner.maxlevel", 15);
    }
    //prm_json.put("tol", 0.001);
    //prm_json.put("maxiter", 30);
    std::function<Vec()> quasi;
    int pidx = 1;
    if (block_size == 2)
	pidx = 0;
    
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);

    Vec x(rhs_loc.size());
    //x=0;
    //Vec crhs(rhs_loc);
    Dune::InverseOperatorResult stat;
    //fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);

    // Initialise the list of parameters with initial, min and max values
    int num_cpr_parameters = 3;
    int num_amg_parameters = 13;
    int num_parameters;
    if (pc_Type == "cpr") {
        num_parameters = num_cpr_parameters + num_amg_parameters;
    } else if (pc_Type == "amg") {
        num_parameters = num_amg_parameters;
    }

    std::string cpr_parameters[][5] = {
        {
            "preconditioner.finesmoother.relaxation",
            "0.9",
            "0",
            "1",
            "double"
        },
        {
            "preconditioner.pre_smooth",
            "1",
            "0",
            "10",
            "int"
        },
        {
            "preconditioner.post_smooth",
            "1",
            "0",
            "10",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.alpha",
            "0.33333333333300003",
            "0",
            "1",
            "double"
        },
        {
            "preconditioner.coarsesolver.preconditioner.relaxation",
            "1",
            "0",
            "1",
            "double"
        },
        {
            "preconditioner.coarsesolver.preconditioner.coarsenTarget",
            "4200",
            "100",
            "10000",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.pre_smooth",
            "1",
            "0",
            "10",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.post_smooth",
            "1",
            "0",
            "10",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.beta",
            "1e-5",
            "0",
            "0.1",
            "double"
        },
        {
            "preconditioner.coarsesolver.preconditioner.skip_isolated", //strongly connected to beta
            "1",
            "0",
            "1",
            "bool"
        },
        {
            "preconditioner.coarsesolver.preconditioner.prolongationdamping",
            "1.6",
            "0.5", //should perhaps be 1
            "3",
            "double"
        },
        {
            "preconditioner.coarsesolver.preconditioner.maxdistance",
            "2",
            "1",
            "10",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.maxconnectivity",
            "15",
            "3",
            "40",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.maxaggsize",
            "6",
            "2",
            "40",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.minaggsize",
            "4",
            "1",
            "20",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.gamma",
            "1",
            "1",
            "3",
            "int"
        }
    };
    std::string amg_parameters[][5] = {
        {
            "preconditioner.alpha",
            "0.33333333333300003",
            "0",
            "1",
            "double"
        },
        {
            "preconditioner.relaxation",
            "1",
            "0",
            "1",
            "double"
        },
        {
            "preconditioner.coarsenTarget",
            "4200",
            "100",
            "10000",
            "int"
        },
        {
            "preconditioner.pre_smooth",
            "1",
            "0",
            "10",
            "int"
        },
        {
            "preconditioner.post_smooth",
            "1",
            "0",
            "10",
            "int"
        },
        {
            "preconditioner.beta",
            "1e-5",
            "0",
            "0.1",
            "double"
        },
        {
            "preconditioner.skip_isolated", //strongly connected to beta
            "1",
            "0",
            "1",
            "bool"
        },
        {
            "preconditioner.prolongationdamping",
            "1.6",
            "0.5", //should perhaps be 1
            "3",
            "double"
        },
        {
            "preconditioner.maxdistance",
            "2",
            "1",
            "10",
            "int"
        },
        {
            "preconditioner.maxconnectivity",
            "15",
            "3",
            "40",
            "int"
        },
        {
            "preconditioner.maxaggsize",
            "6",
            "2",
            "40",
            "int"
        },
        {
            "preconditioner.minaggsize",
            "4",
            "1",
            "20",
            "int"
        },
        {
            "preconditioner.gamma",
            "1",
            "1",
            "3",
            "int"
        }
    };

    std::string preconditioner_parameters[num_parameters][5];
    int counter = 0;
    while (counter < num_parameters) {
        for (int i = 0; i < 5; i++) {
            if (pc_Type == "cpr") {
                preconditioner_parameters[counter][i] = cpr_parameters[counter][i];
            } else if (pc_Type == "amg") {
                preconditioner_parameters[counter][i] = amg_parameters[counter][i];
            }
        }
        counter++;
    }

    // Initialise other necessary variables and arrays
    int num_time_measurement = 3;
    int num_iterations = 30;
    int num_perturbations = 5;
    double times[num_iterations+1];
    double iterations[num_iterations+1];
    double new_parameter_values_list[num_perturbations+1][num_parameters];
    double new_times_list[num_perturbations+1];
    double new_iterations_list[num_parameters+1];
    
    // Set seed for random number generator (do we need this?)
    std::srand(std::time(0));

    // Start main loop over number of iterations (plus one, since we need the initial run)
    for (int i = 0; i < num_iterations + 1; i++) {
        
        // First iterations is just used to solve using the initial parameter values
        if (i == 0) {
            for (int j = 0; j < num_parameters; j++) {
                prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
            }
            auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);
            times[0] = std::numeric_limits<double>::infinity();
            iterations[0] = std::numeric_limits<double>::infinity();
            for (int k = 0; k < num_time_measurement; k++) {
                x = 0;
                Vec crhs(rhs_loc);
                fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);
                if (stat.elapsed < times[0]) {
                    times[0] = stat.elapsed;
                    iterations[0] = stat.iterations;
                }
            }
            if (rank == 0) {
                std::cout << "Initial parameters: (" << times[0] << ", " << iterations[0] << ")\n" << std::endl;
            }
            continue;
        }
        
        // Print out iteration information
        if (rank == 0 && i > 0) {
            std::cout << "\nStarting iteration number: " << i << std::endl;
        }

        // Start inner loop over number of perturbations
        for (int p = 0; p < num_perturbations; p++) {
            if (rank == 0 && i > 0) {
                std::cout << "\tPerturbation number: " << p + 1 << std::endl;
            }

            // Reset preconditioner to current best before perturbing values
            for (int j = 0; j < num_parameters; j++) {
                prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
            }

            // Perturb parameter values and ensure they fall within [min, max]
            for (int j = 0; j < num_parameters; j++) {
                double new_value = prm_json.get<double>(preconditioner_parameters[j][0]);

                // Only change some of the parameter values
                if ((double)std::rand() / RAND_MAX > 0.8) {
                    double min_value = stod(preconditioner_parameters[j][2]);
                    double max_value = stod(preconditioner_parameters[j][3]);
                    if (preconditioner_parameters[j][4] == "bool") {
                        if (new_value == 0) {
                            new_value = 1;
                        }
                        else {
                            new_value = 0;
                        }
                    }
                    else if (new_value == 0) {
                        new_value = (max_value - min_value) * (double)std::rand() / RAND_MAX + min_value;
                    }
                    else if (new_value == 1 && preconditioner_parameters[j][4] == "int") {
                        new_value = (max_value - min_value) * (double)std::rand() / RAND_MAX + min_value;
                    }
                    else {
                        double random_number = (2 * (double)std::rand() / RAND_MAX - 1);
                        new_value += new_value * random_number;
                    }

                    if (new_value < min_value) {
                        new_value = min_value;
                    }
                    if (new_value > max_value) {
                        new_value = max_value;
                    }
                    if (preconditioner_parameters[j][4] == "int") {
                        new_value = std::round(new_value);
                    }
                }
                
                if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.minaggsize" || preconditioner_parameters[j][0] == "preconditioner.minaggsize") {
                    int temp_maxaggsize;
                    if (pc_Type == "cpr") {
                        temp_maxaggsize = prm_json.get<int>("preconditioner.coarsesolver.preconditioner.maxaggsize");
                    } else if (pc_Type == "amg") {
                        temp_maxaggsize = prm_json.get<int>("preconditioner.maxaggsize");
                    }
                    if (new_value > temp_maxaggsize) {
                        new_value = temp_maxaggsize;
                    }
                }
                
                prm_json.put(preconditioner_parameters[j][0], new_value);
                new_parameter_values_list[p][j] = new_value;
            }

            for (int j = 0; j < num_parameters; j++) {
                if (rank == 0) {
                    std::cout << preconditioner_parameters[j][0] << ": " << prm_json.get<std::string>(preconditioner_parameters[j][0]) << std::endl;
                }
            }
            
            // Solve the system using the perturbed parameter values
            try {
                auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);
                new_times_list[p] = std::numeric_limits<double>::infinity();
                new_iterations_list[p] = std::numeric_limits<double>::infinity();
                for (int k = 0; k < num_time_measurement; k++) {
                    x = 0;
                    Vec crhs_new(rhs_loc);
                    fs_json->apply(x, crhs_new, prm_json.get<double>("tol", 0.001), stat);
                    if (stat.elapsed < new_times_list[p]) {
                        if (stat.iterations == 200) {
                            new_times_list[p] = std::numeric_limits<double>::infinity();
                            new_iterations_list[p] = std::numeric_limits<double>::infinity();
                            break;
                        }
                        new_times_list[p] = stat.elapsed;
                        new_iterations_list[p] = stat.iterations;
                    }
                }
            } catch(...) {
                new_times_list[p] = std::numeric_limits<double>::infinity();
                new_iterations_list[p] = std::numeric_limits<double>::infinity();
            }
            
            // Print out time and iteration count for perturbed parameter values
            if (rank == 0) {
                std::cout << "\t(" << new_times_list[p] << ", " << new_iterations_list[p] << ")\n" << std::endl;
            }
        }

        // Print out information about gradient step
        if (rank == 0 && i > 0) {
            std::cout << "\tGradient step" << std::endl;
        }

        // Reset preconditioner to current best before calculating gradient
        for (int j = 0; j < num_parameters; j++) {
            prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
        }

        // Calculate gradient values based on perturbed parameter values and results from
        // solving the linear system with these values
        double array[num_parameters];
        double new_gradient_value;
        for (int j = 0; j < num_parameters; j++) {
            double old_value = stod(preconditioner_parameters[j][1]);
            double min_value = stod(preconditioner_parameters[j][2]);
            double max_value = stod(preconditioner_parameters[j][3]);
            array[j] = 0;
            int num_completed_computations = 0;
            
            for (int p = 0; p < num_perturbations; p++) {
                if (std::isinf(new_times_list[p])) {
                    continue;
                }
                num_completed_computations++;
                array[j] += (times[i-1] - new_times_list[p]) * (new_parameter_values_list[p][j] - old_value);
            }
            //new_gradient_value = old_value + array[j] / (times[i-1] * num_completed_computations);
            if (old_value == 0) {
                new_gradient_value = array[j] / (times[i-1] * num_completed_computations);
            }
            else {
                new_gradient_value = old_value + array[j] / (old_value * times[i-1] * num_completed_computations);
            }
            if (new_gradient_value < min_value) {
                new_gradient_value = min_value;
            }
            if (new_gradient_value > max_value) {
                new_gradient_value = max_value;
            }
            if (preconditioner_parameters[j][4] == "int" || preconditioner_parameters[j][4] == "bool") {
                new_gradient_value = std::round(new_gradient_value);
            }
            if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.minaggsize" || preconditioner_parameters[j][0] == "preconditioner.minaggsize") {
                int temp_maxaggsize;
                if (pc_Type == "cpr") {
                    temp_maxaggsize = prm_json.get<int>("preconditioner.coarsesolver.preconditioner.maxaggsize");
                } else if (pc_Type == "amg") {
                    temp_maxaggsize = prm_json.get<int>("preconditioner.maxaggsize");
                }
                if (new_gradient_value > temp_maxaggsize) {
                    new_gradient_value = temp_maxaggsize;
                }
            }
            new_parameter_values_list[num_perturbations][j] = new_gradient_value;
            prm_json.put(preconditioner_parameters[j][0], new_parameter_values_list[num_perturbations][j]);
        }

        // Ensure that not both pre- and post-smooth are 0
        if (prm_json.get<double>("preconditioner.pre_smooth") == 0 && prm_json.get<double>("preconditioner.post_smooth") == 0) {
            if ((double)std::rand() / RAND_MAX > 0.5) {
                prm_json.put("preconditioner.pre_smooth", 1);
                for (int j = 0; j < num_parameters; j++) {
                    if (preconditioner_parameters[j][0] == "preconditioner.pre_smooth") {
                        new_parameter_values_list[num_perturbations][j] = 1;
                    }
                }
            }
            else {
                prm_json.put("preconditioner.post_smooth", 1);
                for (int j = 0; j < num_parameters; j++) {
                    if (preconditioner_parameters[j][0] == "preconditioner.post_smooth") {
                        new_parameter_values_list[num_perturbations][j] = 1;
                    }
                }
            }
        }

        for (int j = 0; j < num_parameters; j++) {
            if (rank == 0) {
                std::cout << preconditioner_parameters[j][0] << ": " << prm_json.get<std::string>(preconditioner_parameters[j][0]) << std::endl;
            }
        }

        // Solve the system using the gradient values
        try {
            auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);
            new_times_list[num_perturbations] = std::numeric_limits<double>::infinity();
            new_iterations_list[num_perturbations] = std::numeric_limits<double>::infinity();
            for (int k = 0; k < num_time_measurement; k++) {
                x = 0;
                Vec crhs_new(rhs_loc);
                fs_json->apply(x, crhs_new, prm_json.get<double>("tol", 0.001), stat);
                if (stat.elapsed < new_times_list[num_perturbations]) {
                    if (stat.iterations == 200) {
                        new_times_list[num_perturbations] = std::numeric_limits<double>::infinity();
                        new_iterations_list[num_perturbations] = std::numeric_limits<double>::infinity();
                        break;
                    }
                    new_times_list[num_perturbations] = stat.elapsed;
                    new_iterations_list[num_perturbations] = stat.iterations;
                }
            }
        } catch(...) {
                new_times_list[num_perturbations] = std::numeric_limits<double>::infinity();
                new_iterations_list[num_perturbations] = std::numeric_limits<double>::infinity();
        }

        // Print out time and iteration count for gradient parameter values
        if (rank == 0) {
            std::cout << "\t(" << new_times_list[num_perturbations] << ", " << new_iterations_list[num_perturbations] << ")\n" << std::endl;
        }

        // Find the fastest parameter values (from perturbed plus gradient)
        int min_index = 0;
        bool is_gradient_lowest = false;
        for (int indx = 1; indx < num_perturbations + 1; indx++) {
            if (new_times_list[indx] < new_times_list[min_index]) {
                min_index = indx;
            }
        }
        if (min_index == num_perturbations) {
            is_gradient_lowest = true;
        }

        // Compare the fastest parameter values with the current best and update
        // current best if the new one is faster
        if (new_times_list[min_index] < times[i-1]) {
            times[i] = new_times_list[min_index];
            iterations[i] = new_iterations_list[min_index];
            if (rank == 0) {
                std::cout << "Found a better parameter set!" << std::endl;
                std::cout << "Time reduced from " << times[i-1] << " to " << times[i] << std::endl;
                if (is_gradient_lowest) {
                    std::cout << "Gradient update" << std::endl;
                }
                else {
                    std::cout << "Non-gradient update" << std::endl;
                }
            }
            for (int j = 0; j < num_parameters; j++) {
                preconditioner_parameters[j][1] = std::to_string(new_parameter_values_list[min_index][j]);
                if (rank == 0) {
                    std::cout << preconditioner_parameters[j][0] << ": " << preconditioner_parameters[j][1] << std::endl;
                }
            }
        }
        else {
            times[i] = times[i-1];
            iterations[i] = iterations[i-1];
        }
    }

    // Print out final information about times and iterations over the optimisation,
    // and the final best parameter values
    if (rank == 0) {
        std::cout << "\n(Times, iterations):" << std::endl;
        for (int i = 0; i < sizeof(times)/sizeof(times[0]); i++) {
            std::cout << "(" << times[i] << ", " << iterations[i] << ")" << std::endl;
        }
        std::cout << "\nParameter values:" << std::endl;
        for (int i = 0; i < num_parameters; i++) {
            std::cout << preconditioner_parameters[i][0] << ": " << preconditioner_parameters[i][1] << std::endl;
        }
    }
}


template<class Mat, class Vec>
void gen_dim_jsonSolve_extended(int argc, char** argv)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;                 // solveParallel/ghostLastOperations.hpp
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    findZeroDiag(A_loc, rhs_loc);
    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);

    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    prm_json.put("tol", 0.001);
    prm_json.put("maxiter", 30);
    std::function<Vec()> quasi;
    int pidx = 1;

    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);

    Vec x(rhs_loc.size());
    x=0;
    Vec crhs(rhs_loc);
    Dune::InverseOperatorResult stat;
    fs_json->apply(x, crhs, 0.001, stat);

    //flexibleSolverAndPreconTimer(*fs_json, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.0001, true);

    using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

    using Smoother = ILU;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    
    //Criterion criterion(15, 1200,1.2,1.6,static_cast<Dune::Amg::AccumulationMode>(1));//, true);

    Criterion criterion(15, prm_json.get<int>("coarsenTarget", 1200));
    auto pc_child = prm_json.get_child_optional("preconditioner");
    setCrit(criterion, *pc_child );
    criterion.setRandomParallelGhostIndexOrder(false);
    SmootherArgs smootherArgs;
    setOpmILU0args(smootherArgs, *pc_child);
    auto amg = std::make_shared<AMGCPR>(linOp, criterion, smootherArgs, *parComm);

    int verb = 0;
    if (rank == 0) {verb=2;}
    Solver solver(linOp, sp, *amg, 0.01, 30, verb);
    Vec x2(rhs_loc.size());
    Vec crhs2(rhs_loc);
    x2=0;
    Dune::InverseOperatorResult stat2;
    solver.apply(x2, crhs2, stat2);    
}

template<class Mat, class Vec>
void gen_dim_jsonSolve_compare_parallel(int argc, char** argv)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;                 // solveParallel/ghostLastOperations.hpp
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);

    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    prm_json.put("tol", 0.000001);
    prm_json.put("maxiter", 300);
    std::function<Vec()> quasi;
    int pidx = 1;

    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);

    Vec x(rhs_loc.size());
    x=0;
    Vec crhs(rhs_loc);
    Dune::InverseOperatorResult stat;
    fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);

    bool writeSolution = true;
    if (writeSolution) {

	if (cc.size() == 1) {

	    std::ofstream xfile("Solution.vec");
	    xfile.precision (std::numeric_limits<double>::digits10 + 1);
	    Dune::writeMatrixMarket(x, xfile);
	    xfile.close();
	}
	else {
	    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
	    Vec serialSol;
	    readMatMarketObject(serialSol ,"Solution.vec");

	    auto indexSet = parComm->indexSet();

	    double sumDiff = 0;
	    double infNorm = 0;
	    int ownerCells = 0;

	    double pDiff = 0;
	    double satW = 0;
	    double satG = 0;
	    for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {

		if (idx->local().attribute()==AttributeSet::owner) {

		    auto diff = x[idx->local().local()] - serialSol[idx->global()-1];
		    //auto diff = rhs_loc[idx->local().local()] - rhs[idx->global()];
		    
		    sumDiff += diff.two_norm ();
		    //infNorm = std::max(infNorm, diff.infinity_norm ());
		    ownerCells += 1;

		    pDiff += std::abs(diff[1]);
		    satW += std::abs(diff[0]);
		    satG += std::abs(diff[2]);

		    if (idx->global() == 0) {
			std::cout << idx->local().local()<< " " << idx->global()<< " "<< x[idx->local().local()][1]<< " " << serialSol[idx->global()-1][1] << std::endl;
		    }
		    if (false) {
			std::cout << idx->local().local()<< " " << idx->global()<< " "<< x[idx->local().local()][1]<< " " << serialSol[idx->global()-1][1] << std::endl;
		    }
		    
		}
	    }

	    //std::cout<< rank << ": "<< sumDiff<< " " << sumDiff/ownerCells<< " "<< infNorm<< std::endl;

	    std::cout<< rank << ": "<< pDiff<< " " << satW<< " " << satG<< " " << pDiff/ownerCells<< " " << satW/ownerCells<< " " << satG/ownerCells<< std::endl;
	}
    }
}

