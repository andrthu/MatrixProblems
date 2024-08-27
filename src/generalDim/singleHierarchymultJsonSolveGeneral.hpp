/*
  Copyright 2024 Andreas Thune.

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

#ifndef OPM_SINGLEHIERARCHYMULTJSONSOLVEGENERAL_HEADER_INCLUDED
#define OPM_SINGLEHIERARCHYMULTJSONSOLVEGENERAL_HEADER_INCLUDED

#endif // OPM_SINGLEHIERARCHYMULTJSONSOLVEGENERAL_HEADER_INCLUDED


template<class Mat, class Vec>
void gen_dim_jsonSolve_mult_sys_same_hir(std::vector<std::string> systemDirs)
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

    std::vector<Mat> systems;
    std::vector<Vec> rhs;

    std::vector<ScalarProduct> sps;
    
    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    std::vector<int> mpiVec;
    for (int i = 0; i < systemDirs.size(); ++i) {
        if ( boost::algorithm::ends_with( systemDirs[i], ".json") ) {
            DR.dict[12] = systemDirs[i];
        }
        else {
            Mat A_loc;
            Vec rhs_loc;
            mpiVec = readMatOnRootAndDist(systemDirs[i], A_loc, rhs_loc, DR, comm, parComm, cc, mpiVec, true, i!=0);

            systems.push_back(A_loc);
            rhs.push_back(rhs_loc);
            sps.push_back(ScalarProduct(*parComm));
        }
    }

    if (rank == 0) {std::cout << std::endl;}

    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	    prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    if (pc_Type == "cpr") {
        prm_json.put("preconditioner.coarsesolver.preconditioner.maxlevel", 15);
        prm_json.put("preconditioner.coarsesolver.maxiter", 1);
    } else if (pc_Type == "amg") {
        prm_json.put("preconditioner.maxlevel", 15);
    }
    if (pc_Type == "cpr") {
        prm_json.put("maxiter", 30);
    }
    else if (pc_Type == "amg") {
        prm_json.put("maxiter", 150);
    }
    
	Dune::InverseOperatorResult stat;
	
    // Initialise the list of parameters with initial, min and max values
    int num_cpr_parameters = 4; // REMEMBER to change this when adding or removing parameters!
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
            "0.99",
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
            "0",//"1e-5",
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
            "1.0",//"1.6",
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
            "2",
            "20",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.gamma",
            "1",
            "1",
            "2",
            "int"
        },
        {
            "preconditioner.coarsesolver.preconditioner.smoother",
            "0", // REMEMBER TO CHANGE THIS!!!
            "0",
            "4",
            "int"
        }
    };
    std::string amg_parameters[][5] = {
        {
            "preconditioner.alpha",
            "0.33333333333300003",
            "0",
            "0.99",
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
            "0",//"1e-5",
            "0",
            "0.1",
            "double"
        },
        /*{
            "preconditioner.skip_isolated", //strongly connected to beta
            "1",
            "0",
            "1",
            "bool"
        },*/
        {
            "preconditioner.prolongationdamping",
            "1.0",//"1.6",
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
            "2",
            "20",
            "int"
        },
        {
            "preconditioner.gamma",
            "1",
            "1",
            "2",
            "int"
        },
        {
            "preconditioner.smoother",
            "0", // REMEMBER TO CHANGE THIS!!!
            "0",
            "4",
            "int"
        }
    };

    std::string smoothers[] = {
        "ILU0",
        "Jac",
        "GS",
        "SOR",
        "SSOR"
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
    int sysNum = systems.size();
    double max_time_limit = 3;
    int num_time_measurement = 3;
    int num_iterations = 30;
    int num_perturbations = 5;
    double times[num_iterations+1][sysNum][2];
    double iterations[num_iterations+1][sysNum];
    double new_parameter_values_list[num_perturbations+1][num_parameters];
    double new_times_list[num_perturbations+1][sysNum][2];
    double new_iterations_list[num_parameters+1][sysNum];
    double total_best_relative_improvement = 0;
    
    // Set seed for random number generator (do we need this?)
    std::srand(std::time(NULL));
    double random_number;

    // Initialise variables needed to get info about update time (i.e. updating preconditioner before new time step)
    std::ostringstream oss;
    size_t pos = 0;

    std::chrono::steady_clock::time_point start_time_optim = std::chrono::steady_clock::now();

    // Start main loop over number of iterations (plus one, since we need the initial run)
    for (int i = 0; i < num_iterations + 1; i++) {

        // First iterations is just used to solve using the initial parameter values
        if (i == 0) {
            for (int j = 0; j < num_parameters; j++) {
                if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                    preconditioner_parameters[j][0] == "preconditioner.smoother") {
                    prm_json.put(preconditioner_parameters[j][0], smoothers[std::stoi(preconditioner_parameters[j][1])]);
                    continue;
                }
                prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
            }

            if (rank == 0) {
                std::cout << "\nSolving with initial (default) parameters:" << std::endl;
            }

            for (int j = 0; j < num_parameters; j++) {
                if (rank == 0) {
                    std::cout << preconditioner_parameters[j][0] << ": " << prm_json.get<std::string>(preconditioner_parameters[j][0]) << std::endl;
                }
            }

            try {
                Mat* matrix = &systems[0];
                auto glo = std::make_unique<GLO>(*matrix, *parComm);

                // Create QuasiImpesWeights function used for CPR
                int pidx = 1;
                if (block_size == 2) { pidx = 0; }
                std::function<Vec()> quasi;
                quasi = [matrix, pidx]() {
                    return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(*matrix, pidx, false);
                };

                // Create linear solver. AMG hierarchy is set-up here based on the systems[0] matrix 
                cc.barrier();
                auto fs_json = std::make_unique<FlexibleSolverType>(*glo, *parComm, prm_json, quasi, pidx);
                cc.barrier();

                for (int sys = 0; sys < sysNum; sys++) {

                    if (rank == 0) {
                        std::cout << "\n\tSolving linear system number " << sys + 1 << " out of " << sysNum << " systems." << std::endl;
                    }

                    iterations[0][sys] = 0;
                    times[0][sys][0] = std::numeric_limits<double>::infinity();
                    times[0][sys][1] = std::numeric_limits<double>::infinity();

                    matrix = &systems[sys];
                    
                    for (int k = 0; k < num_time_measurement; k++) {
                        double temp_time = 0;
                        double temp_update_time;
                        int temp_iterations;
                        
                        Vec crhs(rhs[sys]);
                        Vec x(crhs.size());
                        x = 0;

                        auto cout_buff = std::cout.rdbuf(oss.rdbuf());

                        cc.barrier();
                        fs_json->preconditioner().update();
                        cc.barrier();
                        
                        std::cout.rdbuf(cout_buff);
                        if (rank == 0) {
                            std::string output = oss.str();
                            pos = output.find("levels ");
                            output.erase(0, pos + 7);
                            pos = output.find(" seconds");
                            output.erase(pos, 9);
                            temp_time += stod(output);
                            temp_update_time = stod(output);
                            std::cout << "\nUpdate time for preconditioner was: " << output << std::endl;
                            oss.str("");
                            oss.clear();
                        }

                        cc.broadcast(&temp_time, 1, 0);
                        cc.broadcast(&temp_update_time, 1, 0);

                        cc.barrier();
                        fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);
                        cc.barrier();

                        temp_time += stat.elapsed;
                        temp_iterations = stat.iterations;

                        cc.broadcast(&temp_time, 1, 0);
                        cc.broadcast(&temp_iterations, 1, 0);

                        if (temp_time < times[0][sys][0]) {
                            if (temp_iterations == prm_json.get<int>("maxiter")) {
                                times[0][sys][0] = std::numeric_limits<double>::infinity();
                                times[0][sys][1] = std::numeric_limits<double>::infinity();
                                iterations[0][sys] = std::numeric_limits<double>::infinity();
                                break;
                            }
                            times[0][sys][0] = temp_time;
                            times[0][sys][1] = temp_update_time;
                            if (k == 0) {
                                iterations[0][sys] = temp_iterations;
                            }
                        }
                    }
                }
            } catch(...) {
                std::cout << "Couldn't solve the linear system with default parameters." << std::endl;
            }
            if (rank == 0) {
                std::cout << "\nInitial times: " << times[0][0][0];
                for (int sys = 1; sys < sysNum; sys++) {
                    std::cout << ", " << times[0][sys][0];
                }
                std::cout << "\nInitial iteration counts: " << iterations[0][0];
                for (int sys = 1; sys < sysNum; sys++) {
                    std::cout << ", " << iterations[0][sys];
                }
                std::cout << "\n\n";
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
                if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                    preconditioner_parameters[j][0] == "preconditioner.smoother") {
                    prm_json.put(preconditioner_parameters[j][0], smoothers[std::stoi(preconditioner_parameters[j][1])]);
                    continue;
                }
                prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
            }

            // Perturb parameter values and ensure they fall within [min, max]
            for (int j = 0; j < num_parameters; j++) {
                double new_value = stod(preconditioner_parameters[j][1]);

                // Only change some of the parameter values
                random_number = (double)std::rand();
                cc.broadcast(&random_number, 1, 0);
                if (random_number / RAND_MAX < 0.6) {
                    double min_value = stod(preconditioner_parameters[j][2]);
                    double max_value = stod(preconditioner_parameters[j][3]);
                    
                    double amount_of_perturbation = std::max(0.3 * new_value, 3.0);
                    double lower_bound = std::max(new_value - amount_of_perturbation, min_value);
                    double upper_bound = std::min(new_value + amount_of_perturbation, max_value);

                    if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.minaggsize" || preconditioner_parameters[j][0] == "preconditioner.minaggsize") {
                        double temp_maxaggsize;
                        if (pc_Type == "cpr") {
                            temp_maxaggsize = prm_json.get<int>("preconditioner.coarsesolver.preconditioner.maxaggsize");
                        } else if (pc_Type == "amg") {
                            temp_maxaggsize = prm_json.get<int>("preconditioner.maxaggsize");
                        }
                        upper_bound = std::min(upper_bound, temp_maxaggsize);
                    }
                    if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.maxaggsize" || preconditioner_parameters[j][0] == "preconditioner.maxaggsize") {
                        double temp_minaggsize;
                        if (pc_Type == "cpr") {
                            temp_minaggsize = prm_json.get<int>("preconditioner.coarsesolver.preconditioner.minaggsize");
                        } else if (pc_Type == "amg") {
                            temp_minaggsize = prm_json.get<int>("preconditioner.minaggsize");
                        }
                        lower_bound = std::max(lower_bound, temp_minaggsize);
                    }

                    random_number = (double)std::rand();
                    cc.broadcast(&random_number, 1, 0);
                    new_value = (upper_bound - lower_bound) * random_number / RAND_MAX + lower_bound;
                    
                    if (preconditioner_parameters[j][4] == "bool") {
                        if (new_value == 0) {
                            new_value = 1;
                        }
                        else {
                            new_value = 0;
                        }
                    }

                    if (preconditioner_parameters[j][4] == "int" || preconditioner_parameters[j][4] == "bool") {
                        new_value = std::round(new_value);
                    }
                }
                                
                if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                    preconditioner_parameters[j][0] == "preconditioner.smoother") {
                    new_parameter_values_list[p][j] = stod(preconditioner_parameters[j][1]);
                    //new_parameter_values_list[p][j] = new_value;
                    //prm_json.put(preconditioner_parameters[j][0], smoothers[int(new_value)]);
                    continue;
                }
                new_parameter_values_list[p][j] = new_value;
                prm_json.put(preconditioner_parameters[j][0], new_value);
            }

            for (int j = 0; j < num_parameters; j++) {
                if (rank == 0) {
                    std::cout << preconditioner_parameters[j][0] << ": " << prm_json.get<std::string>(preconditioner_parameters[j][0]) << std::endl;
                }
            }
            
            // Solve the system using the perturbed parameter values
            try {
                Mat* matrix = &systems[0];
                auto glo = std::make_unique<GLO>(*matrix, *parComm);

                // Create QuasiImpesWeights function used for CPR
                int pidx = 1;
                if (block_size == 2) { pidx = 0; }
                std::function<Vec()> quasi;
                quasi = [matrix, pidx]() {
                    return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(*matrix, pidx, false);
                };

                // Create linear solver. AMG hierarchy is set-up here based on the systems[0] matrix 
                cc.barrier();
                auto fs_json = std::make_unique<FlexibleSolverType>(*glo, *parComm, prm_json, quasi, pidx);
                cc.barrier();

                bool not_converged = false;
                bool too_slow = false;

                for (int sys = 0; sys < sysNum; sys++) {

                    if (not_converged || too_slow) {
                        for (int system = 0; system < sysNum; system++) {
                            new_times_list[p][system][0] = std::numeric_limits<double>::infinity();
                            new_times_list[p][system][1] = std::numeric_limits<double>::infinity();
                            new_iterations_list[p][system] = std::numeric_limits<double>::infinity();
                        }
                        if (rank == 0) {
                            std::cout << "Too slow or did not converge." << std::endl;
                        }
                        break;
                    }

                    if (rank == 0) {
                        std::cout << "\n\tSolving linear system number: " << sys + 1 << " out of " << sysNum << " systems." << std::endl;
                    }

                    new_iterations_list[p][sys] = 0;
                    new_times_list[p][sys][0] = std::numeric_limits<double>::infinity();
                    new_times_list[p][sys][1] = std::numeric_limits<double>::infinity();

                    matrix = &systems[sys];

                    for (int k = 0; k < num_time_measurement; k++) {
                        double temp_time = 0;
                        double temp_update_time;
                        int temp_iterations;

                        Vec crhs(rhs[sys]);
                        Vec x(crhs.size());
                        x = 0;

                        auto cout_buff = std::cout.rdbuf(oss.rdbuf());

                        cc.barrier();
                        fs_json->preconditioner().update();
                        cc.barrier();

                        std::cout.rdbuf(cout_buff);
                        if (rank == 0) {
                            std::string output = oss.str();
                            pos = output.find("levels ");
                            output.erase(0, pos + 7);
                            pos = output.find(" seconds");
                            output.erase(pos, 9);
                            temp_time += stod(output);
                            temp_update_time = stod(output);
                            std::cout << "\nUpdate time for preconditioner was: " << output << std::endl;
                            oss.str("");
                            oss.clear();
                        }

                        cc.broadcast(&temp_time, 1, 0);
                        cc.broadcast(&temp_update_time, 1, 0);

                        cc.barrier();
                        fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);
                        cc.barrier();

                        temp_time += stat.elapsed;
                        temp_iterations = stat.iterations;

                        cc.broadcast(&temp_time, 1, 0);
                        cc.broadcast(&temp_iterations, 1, 0);

                        if (temp_time > max_time_limit * times[0][sys][0]) {
                            too_slow = true;
                            break;
                        }

                        if (temp_iterations == prm_json.get<int>("maxiter")) {
                            not_converged = true;
                            break;
                        }

                        if (temp_time < new_times_list[p][sys][0]) {
                            new_times_list[p][sys][0] = temp_time;
                            new_times_list[p][sys][1] = temp_update_time;
                            if (k == 0) {
                                new_iterations_list[p][sys] = temp_iterations;
                            }
                        }
                    }
                }
            } catch(...) {
                for (int sys = 0; sys < sysNum; sys++) {
                    new_times_list[p][sys][0] = std::numeric_limits<double>::infinity();
                    new_times_list[p][sys][1] = std::numeric_limits<double>::infinity();
                    new_iterations_list[p][sys] = std::numeric_limits<double>::infinity();
                }
            }
            // Print out time and iteration count for perturbed parameter values
            if (rank == 0) {
                std::cout << "\nTimes: " << new_times_list[p][0][0];
                for (int sys = 1; sys < sysNum; sys++) {
                    std::cout << ", " << new_times_list[p][sys][0];
                }
                std::cout << "\nIterations: " << new_iterations_list[p][0];
                for (int sys = 1; sys < sysNum; sys++) {
                    std::cout << ", " << new_iterations_list[p][sys];
                }
                std::cout << "\n\n";
            }
        }

        // Print out information about gradient step
        if (rank == 0 && i > 0) {
            std::cout << "\tGradient step" << std::endl;
        }

        // Reset preconditioner to current best before calculating gradient
        for (int j = 0; j < num_parameters; j++) {
            if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                preconditioner_parameters[j][0] == "preconditioner.smoother") {
                prm_json.put(preconditioner_parameters[j][0], smoothers[std::stoi(preconditioner_parameters[j][1])]);
                continue;
            }
            prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
        }

        // Calculate gradient values based on perturbed parameter values and results from
        // solving the linear system with these values
        double array[num_parameters];
        double new_gradient_value;
        bool all_infinite = false;
        for (int j = 0; j < num_parameters; j++) {
            double old_value = stod(preconditioner_parameters[j][1]);
            double min_value = stod(preconditioner_parameters[j][2]);
            double max_value = stod(preconditioner_parameters[j][3]);
            array[j] = 0;
            int num_completed_computations = 0;
            
            for (int p = 0; p < num_perturbations; p++) {
                bool is_infinite = false;
                for (int sys = 0; sys < sysNum; sys++) {
                    if (std::isinf(new_times_list[p][sys][0])) {
                        is_infinite = true;
                        break;
                    }
                }
                if (is_infinite) {
                    continue;
                }
                num_completed_computations++;
                for (int sys = 0; sys < sysNum; sys++) {
                    array[j] += (times[i-1][sys][0] - new_times_list[p][sys][0]) * (new_parameter_values_list[p][j] - old_value);
                }
            }

            if (num_completed_computations == 0) {
                all_infinite = true;
                break;
            }

            double total_time = 0;
            for (int sys = 0; sys < sysNum; sys++) {
                total_time += times[i-1][sys][0];
            }

            //new_gradient_value = old_value + array[j] / (total_time * num_completed_computations);
            if (old_value == 0) {
                new_gradient_value = array[j] / (total_time * num_completed_computations);
            }
            else {
                new_gradient_value = old_value + array[j] / (old_value * total_time * num_completed_computations);
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
            if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                preconditioner_parameters[j][0] == "preconditioner.smoother") {
                prm_json.put(preconditioner_parameters[j][0], smoothers[int(new_gradient_value)]);
                continue;
            }
            prm_json.put(preconditioner_parameters[j][0], new_gradient_value);
        }

        if (all_infinite) {
            i--;
            continue;
        }

        // Ensure that not both pre- and post-smooth are 0
        if (prm_json.get<double>("preconditioner.pre_smooth") == 0 && prm_json.get<double>("preconditioner.post_smooth") == 0) {
            random_number = (double)std::rand();
            cc.broadcast(&random_number, 1, 0);
            if (random_number / RAND_MAX > 0.5) {
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
        if (pc_Type == "cpr") {
            if (prm_json.get<double>("preconditioner.coarsesolver.preconditioner.pre_smooth") == 0 &&
                prm_json.get<double>("preconditioner.coarsesolver.preconditioner.post_smooth") == 0) {
                random_number = (double)std::rand();
                cc.broadcast(&random_number, 1, 0);
                if (random_number / RAND_MAX > 0.5) {
                    prm_json.put("preconditioner.coarsesolver.preconditioner.pre_smooth", 1);
                    for (int j = 0; j < num_parameters; j++) {
                        if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.pre_smooth") {
                            new_parameter_values_list[num_perturbations][j] = 1;
                        }
                    }
                }
                else {
                    prm_json.put("preconditioner.coarsesolver.preconditioner.post_smooth", 1);
                    for (int j = 0; j < num_parameters; j++) {
                        if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.post_smooth") {
                            new_parameter_values_list[num_perturbations][j] = 1;
                        }
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
            Mat* matrix = &systems[0];
            auto glo = std::make_unique<GLO>(*matrix, *parComm);

            // Create QuasiImpesWeights function used for CPR
            int pidx = 1;
            if (block_size == 2) { pidx = 0; }
            std::function<Vec()> quasi;
            quasi = [matrix, pidx]() {
                return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(*matrix, pidx, false);
            };

            // Create linear solver. AMG hierarchy is set-up here based on the systems[0] matrix 
            cc.barrier();
            auto fs_json = std::make_unique<FlexibleSolverType>(*glo, *parComm, prm_json, quasi, pidx);
            cc.barrier();

            bool not_converged = false;
            bool too_slow = false;

            for (int sys = 0; sys < sysNum; sys++) {
                
                if (not_converged || too_slow) {
                    for (int system = 0; system < sysNum; system++) {
                        new_times_list[num_perturbations][system][0] = std::numeric_limits<double>::infinity();
                        new_times_list[num_perturbations][system][1] = std::numeric_limits<double>::infinity();
                        new_iterations_list[num_perturbations][system] = std::numeric_limits<double>::infinity();
                    }
                    if (rank == 0) {
                        std::cout << "Too slow or did not converge." << std::endl;
                    }
                    break;
                }
                
                if (rank == 0) {
                    std::cout << "\n\tSolving linear system number: " << sys + 1 << " out of " << sysNum << " systems." << std::endl;
                }

                new_iterations_list[num_perturbations][sys] = 0;
                new_times_list[num_perturbations][sys][0] = std::numeric_limits<double>::infinity();
                new_times_list[num_perturbations][sys][1] = std::numeric_limits<double>::infinity();

                matrix = &systems[sys];

                for (int k = 0; k < num_time_measurement; k++) {
                    double temp_time = 0;
                    double temp_update_time;
                    int temp_iterations;

                    Vec crhs(rhs[sys]);
                    Vec x(crhs.size());
                    x = 0;

                    auto cout_buff = std::cout.rdbuf(oss.rdbuf());

                    cc.barrier();
                    fs_json->preconditioner().update();
                    cc.barrier();

                    std::cout.rdbuf(cout_buff);
                    if (rank == 0) {
                        std::string output = oss.str();
                        pos = output.find("levels ");
                        output.erase(0, pos + 7);
                        pos = output.find(" seconds");
                        output.erase(pos, 9);
                        temp_time += stod(output);
                        temp_update_time = stod(output);
                        std::cout << "\nUpdate time for preconditioner was: " << output << std::endl;
                        oss.str("");
                        oss.clear();
                    }

                    cc.broadcast(&temp_time, 1, 0);
                    cc.broadcast(&temp_update_time, 1, 0);

                    cc.barrier();
                    fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);
                    cc.barrier();

                    temp_time += stat.elapsed;
                    temp_iterations = stat.iterations;

                    cc.broadcast(&temp_time, 1, 0);
                    cc.broadcast(&temp_iterations, 1, 0);

                    if (temp_time > max_time_limit * times[0][sys][0]) {
                        too_slow = true;
                        break;
                    }

                    if (temp_iterations == prm_json.get<int>("maxiter")) {
                        not_converged = true;
                        break;
                    }

                    if (temp_time < new_times_list[num_perturbations][sys][0]) {
                        new_times_list[num_perturbations][sys][0] = temp_time;
                        new_times_list[num_perturbations][sys][1] = temp_update_time;
                        if (k == 0) {
                            new_iterations_list[num_perturbations][sys] = temp_iterations;
                        }
                    }
                }
            }
        } catch(...) {
            for (int sys = 0; sys < sysNum; sys++) {
                new_times_list[num_perturbations][sys][0] = std::numeric_limits<double>::infinity();
                new_times_list[num_perturbations][sys][1] = std::numeric_limits<double>::infinity();
                new_iterations_list[num_perturbations][sys] = std::numeric_limits<double>::infinity();
            }
        }

        // Print out time and iteration count for gradient parameter values
        if (rank == 0) {
            std::cout << "\nGradient times: " << new_times_list[num_perturbations][0][0];
            for (int sys = 1; sys < sysNum; sys++) {
                std::cout << ", " << new_times_list[num_perturbations][sys][0];
            }
            std::cout << "\nGradient iterations: " << new_iterations_list[num_perturbations][0];
            for (int sys = 1; sys < sysNum; sys++) {
                std::cout << ", " << new_iterations_list[num_perturbations][sys];
            }
            std::cout << "\n\n";
        }

        // Find the fastest parameter values (from perturbed plus gradient)
        double best_relative_improvement = -std::numeric_limits<double>::infinity();
        int min_index;
        bool is_gradient_lowest = false;
        for (int indx = 1; indx < num_perturbations + 1; indx++) {
            bool is_infinite = false;
            double relative_improvement = 0;
            for (int sys = 0; sys < sysNum; sys++) {
                if (std::isinf(new_times_list[indx][sys][0])) {
                    is_infinite = true;
                    break;
                }
                relative_improvement += (times[0][sys][0] - new_times_list[indx][sys][0]) / times[0][sys][0];
            }
            if (is_infinite) {
                continue;
            }
            if (relative_improvement > best_relative_improvement) {
                best_relative_improvement = relative_improvement;
                min_index = indx;
            }
        }
        if (min_index == num_perturbations) {
            is_gradient_lowest = true;
        }

        // Compare the fastest parameter values with the current best and update
        // current best if the new one is faster
        if (best_relative_improvement > total_best_relative_improvement) {
            total_best_relative_improvement = best_relative_improvement;
            for (int sys = 0; sys < sysNum; sys++) {
                times[i][sys][0] = new_times_list[min_index][sys][0];
                times[i][sys][1] = new_times_list[min_index][sys][1];
                iterations[i][sys] = new_iterations_list[min_index][sys];
            }

            double total_time[] = { 0, 0 };
            double total_solver_time[] = { 0, 0 };
            double total_update_time[] = { 0, 0 };
            int total_iterations[] = { 0, 0 };
            for (int sys = 0; sys < sysNum; sys++) {
                total_time[0] += times[i-1][sys][0];
                total_time[1] += times[i][sys][0];
                total_solver_time[0] += times[i-1][sys][0] - times[i-1][sys][1];
                total_solver_time[1] += times[i][sys][0] - times[i][sys][1];
                total_update_time[0] += times[i-1][sys][1];
                total_update_time[1] += times[i][sys][1];
                total_iterations[0] += iterations[i-1][sys];
                total_iterations[1] += iterations[i][sys];
            }

            if (rank == 0) {
                std::cout << "\nFound a better parameter set!" << std::endl;
                std::cout << "Time changed from " << total_time[0] << " to " << total_time[1] << std::endl;
                std::cout << "Solver time changed from " << total_solver_time[0] << " to " << total_solver_time[1] << std::endl;
                std::cout << "Update time changed from " << total_update_time[0] << " to " << total_update_time[1] << std::endl;
                std::cout << "Iteration count changed from " << total_iterations[0] << " to " << total_iterations[1] << std::endl;
                if (is_gradient_lowest) {
                    std::cout << "Gradient update" << std::endl;
                }
                else {
                    std::cout << "Non-gradient update" << std::endl;
                }
                std::cout << "Single system times: " << times[i][0][0];
                for (int sys = 1; sys < sysNum; sys++) {
                    std::cout << ", " << times[i][sys][0];
                }
                std::cout << "\nSingle system iterations: " << iterations[i][0];
                for (int sys = 1; sys < sysNum; sys++) {
                    std::cout << ", " << iterations[i][sys];
                }
                std::cout << std::endl;
            }
            for (int j = 0; j < num_parameters; j++) {
                if (preconditioner_parameters[j][4] == "int" || preconditioner_parameters[j][4] == "bool") {
                    preconditioner_parameters[j][1] = std::to_string((int)new_parameter_values_list[min_index][j]);
                }
                else {
                    preconditioner_parameters[j][1] = std::to_string(new_parameter_values_list[min_index][j]);
                }
                if (preconditioner_parameters[j][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                    preconditioner_parameters[j][0] == "preconditioner.smoother") {
                    prm_json.put(preconditioner_parameters[j][0], smoothers[std::stoi(preconditioner_parameters[j][1])]);
                }
                else {
                    prm_json.put(preconditioner_parameters[j][0], preconditioner_parameters[j][1]);
                }
                if (rank == 0) {
                    std::cout << preconditioner_parameters[j][0] << ": " << prm_json.get<std::string>(preconditioner_parameters[j][0]) << std::endl;
                }
            }
        }
        else {
            for (int sys = 0; sys < sysNum; sys++) {
                times[i][sys][0] = times[i-1][sys][0];
                times[i][sys][1] = times[i-1][sys][1];
                iterations[i][sys] = iterations[i-1][sys];
            }
        }
    }

    // Print out final information about times and iterations over the optimisation,
    // and the final best parameter values

    std::chrono::steady_clock::time_point end_time_optim = std::chrono::steady_clock::now();
    double total_time;
    int total_iterations;
    if (rank == 0) {
        std::cout << "\n(Times, iterations):" << std::endl;
        for (int i = 0; i < num_iterations + 1; i++) {
            total_time = 0;
            total_iterations = 0;
            for (int sys = 0; sys < sysNum; sys++) {
                total_time += times[i][sys][0];
                total_iterations += iterations[i][sys];
            }
            std::cout << "(" << total_time << ", " << total_iterations << ")" << std::endl;
        }
        std::cout << "\nParameter values:" << std::endl;
        for (int i = 0; i < num_parameters; i++) {
            if (preconditioner_parameters[i][0] == "preconditioner.coarsesolver.preconditioner.smoother" ||
                preconditioner_parameters[i][0] == "preconditioner.smoother") {
                prm_json.put(preconditioner_parameters[i][0], smoothers[std::stoi(preconditioner_parameters[i][1])]);
            }
            else {
                prm_json.put(preconditioner_parameters[i][0], preconditioner_parameters[i][1]);
            }
            std::cout << preconditioner_parameters[i][0] << ": " << prm_json.get<std::string>(preconditioner_parameters[i][0]) << std::endl;
        }
        std::cout << "Time spent on the optimisation: " <<
                        std::chrono::duration_cast<std::chrono::seconds>(end_time_optim - start_time_optim).count() << " seconds" << std::endl;
    }
}
