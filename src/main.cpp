#include "input_reader.h"
#include "root_function.h"
#include "output_writer.h"

typedef struct{
    DM da;
    std::unique_ptr<RootFunction> rf_wrap;
}   AppCtx;

//?? calculation with 1 core
int main(int argc, char *argv[])
{
    //Initialize MPI
    char help[0];
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));

    //Retrieve information about process and communicator
	int rank;
	int comm_size;
	PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
	PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&comm_size));

    //Create world object for broadcasting complex data types
    boost::mpi::communicator world;
    AppCtx Ctx;

    //Declare objects to be defined by process 0
	std::shared_ptr<CellConfiguration> cc = std::make_shared<CellConfiguration>();
	std::shared_ptr<InputReader> ir = std::make_shared<InputReader>();
    Ctx.rf_wrap = std::make_unique<RootFunction>(cc);
    int nodes_per_process[comm_size]; //Distribution of solution vector nodes over the processes, per unit 
    int solution_vector_length;
    std::string input_folder;
    std::string output_folder;
    PetscInt number_of_variables;
    std::vector<int> node_layer_index;
    
    if(rank == 0)
    {
        if(argc == 1){
            std::cerr << "Error: No input provided when executing nernst" << std::endl;
        }

        //Parse command-line arguments
        std::ostringstream strs;
        std::string parse;

        int arg_i = 1; //?? check if both input arguments provided
        while(arg_i < argc){
            strs << argv[arg_i];
            parse = strs.str();

            //Input folder
            if(parse == "-i"){
                arg_i++;
                if(arg_i < argc){
                    strs.str(std::string());
                    strs << argv[arg_i];
                    input_folder = strs.str();
                }
                else{
                    break;
                }
            }

            //Output folder
            if(parse == "-o"){
                arg_i++;
                if(arg_i < argc){
                    strs.str(std::string());
                    strs << argv[arg_i];
                    output_folder = strs.str();
                }
                else{
                    break;
                }
            }

            strs.str(std::string());
            arg_i++;
        }

        std::string input_file = input_folder + "input.txt"; //?? file existence

        //Create cell configuration object
        cc->construct_class();

        //Create input reader object, add the layers to the cell configuration object based on the input_file
        ir->construct_class(input_file,cc);

        //Show cell configuration to user
        cc->print_configuration();
    }
    //Broadcast information from process 0 to all processes
    boost::mpi::broadcast(world,*cc,0);  //Boost::mpi for complex data types (requires ~1e-4 s)
    boost::mpi::broadcast(world,*ir,0);
    // boost::mpi::broadcast(world,*Ctx.rf_wrap,0);   

    //Distribute solution vector over processes
    cc->calculate_solution_vector_unit_lengths(node_layer_index,solution_vector_length, number_of_variables);
        
    for(int j = 0; j < comm_size; j++){
        nodes_per_process[j] = std::floor(solution_vector_length/comm_size);
    }
    nodes_per_process[comm_size-1] += solution_vector_length - std::floor(solution_vector_length/comm_size)*comm_size;
    
    MPI_Bcast(&solution_vector_length,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(nodes_per_process,comm_size,MPI_INT,0,PETSC_COMM_WORLD); //MPI for simple data types
    MPI_Bcast(&number_of_variables,1,MPI_INT,0,PETSC_COMM_WORLD); //MPI for simple data types
    // MPI_Bcast(&nodes_pp,comm_size,MPI_INT,0,PETSC_COMM_WORLD); //MPI for simple data types

    // Create data management object for distributed array, based on number of nodes per process
    // DM da;
    PetscCall(DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, solution_vector_length, number_of_variables, 1, nodes_per_process, &Ctx.da));    
    PetscCall(DMSetFromOptions(Ctx.da));
    PetscCall(DMSetUp(Ctx.da));

    Ctx.rf_wrap->construct_class(ir->get_double_input("T"), ir->get_double_input("area"),ir->get_bool_input("galvanostatic"), 
                        ir->get_double_input("applied_I_or_V"),
                        comm_size, solution_vector_length, nodes_per_process, node_layer_index);


    // //Create vectors from DMDA
    Vec x;
    PetscCall(DMCreateGlobalVector(Ctx.da, &x)); //Solution vector

    Ctx.rf_wrap->construct_class_per_process(rank,x,Ctx.da);

    // //Set initial guess
    PetscCall(Ctx.rf_wrap->set_initial_guess(x,Ctx.da));

    //Create non-linear solver context
    SNES snes;
    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));

    // //Set root function
    PetscCall(SNESSetDM(snes, Ctx.da));
    PetscCall(SNESSetFunction(snes, NULL, root_function_wrapper, &Ctx));

    // //Set Jacobian
    Mat J;
    PetscCall(MatCreate(PETSC_COMM_WORLD,&J));
    PetscCall(MatSetType(J,MATMPIAIJ));
    PetscCall(MatSetSizes(J, number_of_variables*nodes_per_process[rank],number_of_variables*nodes_per_process[rank],
                          solution_vector_length*number_of_variables,solution_vector_length*number_of_variables));
    PetscCall(MatMPIAIJSetPreallocation(J,solution_vector_length*number_of_variables,NULL,
                                        solution_vector_length*number_of_variables,NULL)); //?? improve this!
    PetscCall(SNESSetJacobian(snes, J, J, jacobian_wrapper, &Ctx));
    PetscCall(SNESSetTolerances(snes, 1e-20, 1e-16, 1e-12, 100, -1));

    PetscCall(SNESSetFromOptions(snes));
    //Solve non-linear system
    auto start = std::chrono::system_clock::now();

    PetscCall(SNESSolve(snes, NULL, x));
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time = end - start;
    std::cout << "Completed in " << (boost::format("%8.6f s ") % elapsed_time.count()).str() << std::endl;
    






    //Check if solver is converged
    if(rank == 0){
        SNESConvergedReason snes_convergence;
        PetscCall(SNESGetConvergedReason(snes,&snes_convergence));
        if(snes_convergence == 0){
            std::cout << "Solver still iterating" << std::endl;
        }
        else if(snes_convergence > 0){
            std::cout << "Solver converged" << std::endl;
        }
        else{
            std::cerr << "Warning: solver not converged (reason: " << Ctx.rf_wrap->get_petsc_convergence_reasons(snes_convergence) << ")" << std::endl;
        }
    }

    //?? out object?
    OutputWriter out(rank,output_folder,input_folder);
    out.write_solution_to_csv(  nodes_per_process,comm_size,solution_vector_length,
                                cc,x,Ctx.rf_wrap->Fluxes,Ctx.rf_wrap->Source_term,Ctx.da);

    PetscCall(SNESDestroy(&snes));
    PetscCall(DMDestroy(&Ctx.da));
    PetscCall(VecDestroy(&x));
    PetscCall(MatDestroy(&J));
    Ctx.rf_wrap.reset(); //Release ownership of unique ptr: call destructor

    PetscCall(PetscFinalize());
    return 0;
}








    ///////////////////////////////////////////////////////////

    // Mat J;
    // PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nodes_per_process[rank], nodes_per_process[rank], solution_vector_length, solution_vector_length, 2, NULL, 2, NULL, &J));
    
    // PetscCall(DMSetMatType(da,MATAIJ));
    // PetscCall(DMCreateMatrix(da,&J));     

    //PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, nodes_per_process[rank], nodes_per_process[rank], solution_vector_length, solution_vector_length, nodes_per_process[rank], NULL, nodes_per_process[rank], NULL, &J));


    //PetscCall(SNESSetJacobian(snes,J,J,SNESComputeJacobianDefaultColor,NULL));
    //PetscCall(SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,NULL));

    
    //PetscCall(SNESComputeJacobian(snes,x,J,NULL));


    ///////////////////////////////////////////////////////////


    // auto start = std::chrono::system_clock::now();


    // auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_time = end - start;
    // std::cout << "Completed in " << (boost::format("%8.6f s ") % elapsed_time.count()).str() << std::endl;




    ///////////////////////////////////////////////////////////


    // PetscInt its;
    // PetscCall(SNESGetIterationNumber(snes, &its));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Number of SNES iterations = %" PetscInt_FMT "\n", its));


    ///////////////////////////////////////////////////////////
