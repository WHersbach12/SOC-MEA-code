#include "root_function.h"


RootFunction::RootFunction(const std::shared_ptr<CellConfiguration>& _cell_configuration):
cell_configuration(_cell_configuration){};

RootFunction::~RootFunction(){
    delete[] this->root_sub_functions;
    delete[] this->layer_index_mapping;

    const int error_code_f = VecDestroy(&f);
    const int error_code_df = VecDestroy(&df);

    if(error_code_f > 0 || error_code_df > 0){
        std::cerr << "Warning: an error was encountered during the destruction of Vec f or Vec df (RootFunction)" << std::endl;
    }
    
}

void RootFunction::construct_class( const double _temperature,
                                    const double _area,
                                    const bool _galvanostatic,
                                    const double _applied_I_or_V,
                                    const int _comm_size,
                                    const int _solution_vector_length,
                                    int* _nodes_per_process,
                                    const std::vector<int> _node_layer_index){
    
    //Write to private fields
    this->temperature = _temperature;
    this->area = _area;
    this->galvanostatic = _galvanostatic;
    this->applied_I_or_V = _applied_I_or_V; 
    this->comm_size = _comm_size;
    this->solution_vector_length = _solution_vector_length;
    this->nodes_per_process = _nodes_per_process;
    this->node_layer_index = _node_layer_index;

    

    //Fill in petsc_convergence_reasons map
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(0,"SNES_CONVERGED_ITERATING"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(2,"SNES_CONVERGED_FNORM_ABS"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(3,"SNES_CONVERGED_FNORM_RELATIVE"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(4,"SNES_CONVERGED_SNORM_RELATIVE"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(5,"SNES_CONVERGED_ITS"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(6,"SNES_BREAKOUT_INNER_ITER"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-1,"SNES_DIVERGED_FUNCTION_DOMAIN"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-2,"SNES_DIVERGED_FUNCTION_COUNT"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-3,"SNES_DIVERGED_LINEAR_SOLVE"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-4,"SNES_DIVERGED_FNORM_NAN"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-5,"SNES_DIVERGED_MAX_IT"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-6,"SNES_DIVERGED_LINE_SEARCH"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-7,"SNES_DIVERGED_INNER"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-8,"SNES_DIVERGED_LOCAL_MIN"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-9,"SNES_DIVERGED_DTOL"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-10,"SNES_DIVERGED_JACOBIAN_DOMAIN"));
    this->petsc_convergence_reasons.insert(std::make_pair<int,std::string>(-11,"SNES_DIVERGED_TR_DELTA"));

};


PetscErrorCode RootFunction::construct_class_per_process(const int _rank, Vec x, DM da){
    
    //Set process specific properties
    this->rank = _rank;

    //Duplicate vectors
    PetscCall(VecDuplicate(x,&this->f));
    PetscCall(VecDuplicate(x,&this->df));
    PetscCall(VecDuplicate(x,&this->Fluxes));
    PetscCall(VecDuplicate(x,&this->Source_term));

    //Construct the factory object (local scope only)
    RootSubFunctionFactory rsf_factory;

    //Determine which layers are present for the current process
    const int layer_count = this->cell_configuration->get_layer_count();
    bool layer_present[layer_count]; //Indicates for every layer whether it is present in the current process
    for(int i = 0; i < layer_count; i++){layer_present[i] = 0;} //Set all values to zero

    //Iterate over all units the process is responsible for
    int Process_start_index[this->comm_size];
    Process_start_index[0] = 0;
    for(int i = 1; i < this->comm_size; i++){
        Process_start_index[i] = Process_start_index[i-1] + this->nodes_per_process[i-1];
    }
    for(int i = Process_start_index[this->rank]; i < Process_start_index[this->rank] + this->nodes_per_process[this->rank]; i++){ 
        const int layer_index = this->node_layer_index[i];

        layer_present[layer_index] = 1;

        if(i == Process_start_index[this->rank] && i > 0){ //For the first unit of the process, check if the preceding unit is in the same layer
            const int preceding_layer_index = this->node_layer_index[i-1];
            if(layer_index != preceding_layer_index){ //If the preceding unit is in a different layer, include the preceding layer in the layers present
                                                      //(for the boundary conditions)
                layer_present[preceding_layer_index] = 1;
            }
        }
    }
    //Set the size for the root_sub_functions array
    this->root_sub_functions = new std::shared_ptr<RootSubFunction>[layer_count];

    bool After_dense = 0;
    for(int i = 0; i < layer_count; i++){ //Iterate over all layers
        this->root_sub_functions[i] = rsf_factory.factory_method(   this->rank, this->temperature,
                                                                    this->cell_configuration->get_layer_name(i),
                                                                    this->cell_configuration->get_layer_x_center(i),
                                                                    this->cell_configuration->get_layer_x_face(i),
                                                                    this->cell_configuration->get_layer_gas_species(i),
                                                                    this->galvanostatic, this->applied_I_or_V, this->area, 
                                                                    i, this->cell_configuration);
        if(i > 0){ //If the process is responsible for multiple layers (or the interface between two layers), give layer j access to the object of layer j-1 
        //(for the boundary conditions)
            this->root_sub_functions[i]->set_preceding_layer(this->root_sub_functions[i-1]);
        }
        if(this->cell_configuration->get_layer_is_dense(i)){
            After_dense = 1;
        }
        else if(After_dense){//Last layer - set pointer to fuel gas_compartment
            this->root_sub_functions[i]->set_gas_compartment(this->cell_configuration->get_gas_compartment(0));
        }   
        else{ //First layer - set pointer to air gas_compartment
            this->root_sub_functions[i]->set_gas_compartment(this->cell_configuration->get_gas_compartment(1));
        }

    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode RootFunction::set_initial_guess(Vec x, DM da) const{   
    PetscScalar**   xx;
    PetscInt        xx_start, xx_length, number_of_variables;
    PetscCall(DMDAGetInfo(da,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,&number_of_variables,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE));

    PetscCall(DMDAVecGetArrayDOF(da,x,&xx));

    PetscCall(DMDAGetCorners(da,&xx_start, NULL, NULL, &xx_length, NULL, NULL));

    //Initialize the vector by setting everything to 0.
    for(int i = xx_start; i < xx_start + xx_length; i++){
        for(int j =0; j < number_of_variables; j++){xx[i][j] = 0.0;}
    }
    

    for(int i = xx_start; i < xx_start + xx_length; i++){
        const int layer_index = this->node_layer_index[i];
        const bool is_dense = this->cell_configuration->get_layer_is_dense(layer_index);
        const std::pair<int,int> charge_carriers = this->cell_configuration->get_layer_charge_carriers(layer_index);
        const int number_of_gas_species = this->cell_configuration->get_layer_number_of_gas_species(layer_index);
        const bool is_air_electrode = this->cell_configuration->get_layer_is_air_electrode(layer_index);
        // const std::vector<double> inlet_mole_fractions;

        //Gasses
        for(int j = 0; j < number_of_gas_species; j++){ //Change number of charge to be appropriate for each layer
            // if(is_air_electrode){
                const std::vector<double> inlet_mole_fractions = this->cell_configuration->get_gas_compartment_inlet_mole_fractions(is_air_electrode);
            // }
            // else{
                // const std::vector<double> inlet_mole_fractions = this->cell_configuration->get_gas_compartment_inlet_mole_fractions(!is_air_electrode);
            // }
            const double pressure = this->cell_configuration->get_gas_compartment_inlet_pressure(is_air_electrode);
            const double concentration = pressure/(8.3145*this->temperature);    

            xx[i][j] = concentration*inlet_mole_fractions[j];
            //?? Not working: Porous porous interfaces
            // Porous - porous interfaces
            if(number_of_gas_species < this->cell_configuration->electron_index){
                for(int j = 0; j < number_of_gas_species; j++){ //Change number of charge to be appropriate for each layer
                    xx[i][number_of_gas_species+j] = concentration*inlet_mole_fractions[j];
                }
            }
        }
        

        //Electrons
        for(int j = 0; j < charge_carriers.first; j++){
            if(is_air_electrode && !this->galvanostatic){
                xx[i][this->cell_configuration->electron_index+j] = this->applied_I_or_V*96485.33;
                //First (dominant) electron conductor && potentiostatic mode && air electrode
            }
            else{
                xx[i][this->cell_configuration->electron_index+j] = 0.01*96485.33;
            }
            //Fuel electrode: initial guess always 0 (electrochemical potential electrons always fixed at 0)
            //Galvanostatic: electrochemical potential not fixed - no value to be used for initial guess
            //Secondary electrons (j > 0): not considered to be connected to current collector

        }

        //Oxygen anions
        for(int j = 0; j < charge_carriers.second; j++){
            xx[i][this->cell_configuration->ion_index + j] = 2*96485.33;
        }
    }

    PetscCall(DMDAVecRestoreArrayDOF(da,x,&xx));
    PetscFunctionReturn(0);    
};


PetscErrorCode RootFunction::calculate_root_function(Vec x, Vec f, DM da) const{ //x is the solution vector, f is the function vector

    //Retrieve local vector
    Vec local_x;
    PetscScalar** xx;
    PetscScalar** ff;
    PetscInt xx_start,xx_length, number_of_variables;

    const int electron_index = this->cell_configuration->electron_index;
    const int ion_index = this->cell_configuration->ion_index;

    // PetscCall(DMDAVecGetArrayDOF(da,x,&xx)); Single core, doesnt work for multi
    PetscCall(DMGetLocalVector(da,&local_x));
    PetscCall(DMGlobalToLocalBegin(da,x,INSERT_VALUES, local_x));
    PetscCall(DMGlobalToLocalEnd(da,x,INSERT_VALUES,local_x));

    PetscCall(DMDAVecGetArrayDOF(da,local_x,&xx));
    PetscCall(DMDAVecGetArrayDOF(da,f,&ff));
    PetscCall(DMDAGetInfo(da,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,&number_of_variables,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE));
    PetscCall(DMDAGetCorners(da,&xx_start, NULL, NULL, &xx_length, NULL, NULL));

    //Set all values of ff to zero (ff[j] = flux[j], ff[j] += source[j] -> ff[j] must be zero if there is no flux)
    for(int i = xx_start; i < xx_start + xx_length; i++){
        for(int j = 0; j < number_of_variables; j++){ff[i][j] = 0.0;} //Set all values to zero
    }
    PetscScalar     Western_flux[xx_length][number_of_variables];
    PetscScalar     Eastern_flux[xx_length][number_of_variables];
    PetscScalar     Source[xx_length][number_of_variables];

    //Iterate over all nodes
    for(int i = 0; i < xx_length; i++){
        const int function_index = node_layer_index[xx_start+i];
        const int layer_start_index = this->cell_configuration->layer_start_index[function_index];
        //Calculate western fluxes
        for(int j = 0; j < number_of_variables; j++){
            Western_flux[i][j] = 0.0;
            Eastern_flux[i][j] = 0.0;
            Source[i][j] = 0.0; // Set all values to zero
        }
        // //The flux is calculated in the positive x direction (J = - D dc/dx)
        if(i == 0 && xx_start == 0){
            this->root_sub_functions[function_index]->calculate_flux(xx[xx_start+i],xx_start + i - layer_start_index,Western_flux[i],electron_index,ion_index);
        }
        else{
            this->root_sub_functions[function_index]->calculate_flux(   xx[xx_start+i],xx[xx_start+i-1],xx_start + i - layer_start_index,
                                                                        Western_flux[i],electron_index,ion_index);
        }
    }

    for(int i = 0; i < xx_length-1; i++){
        for(int j = 0; j < number_of_variables; j++){
            // Set corresponding eastern fluxes
            Eastern_flux[i][j] = Western_flux[i+1][j];
        }
    }

    //Final fluxes
    if(this->rank == this->comm_size-1){
        const int function_index = node_layer_index[solution_vector_length-1];
        const int layer_start_index = this->cell_configuration->layer_start_index[function_index]; 
        this->root_sub_functions[function_index]->calculate_flux(   xx[solution_vector_length-1],solution_vector_length - 1 - layer_start_index,
                                                                    Eastern_flux[xx_length-1], electron_index, ion_index);
    }
    else{
        const int function_index = node_layer_index[xx_start + xx_length];
        const int layer_start_index = this->cell_configuration->layer_start_index[function_index];
        this->root_sub_functions[function_index]->calculate_flux(   xx[xx_start + xx_length],xx[xx_start + xx_length-1],xx_start + xx_length - layer_start_index,
                                                                    Eastern_flux[xx_length-1],electron_index,ion_index);
    }

    for(int i = 0; i < xx_length; i++){
        const int function_index = node_layer_index[xx_start+ i];
        this->root_sub_functions[function_index]->calculate_source_term(xx[xx_start + i], Source[i], electron_index, ion_index);
    }

    for(int i = 0; i < xx_length; i++){
        const int function_index = node_layer_index[xx_start + i];
        std::vector<double> x_face = this->cell_configuration->get_layer_x_face(function_index);
        const int layer_start_index = this->cell_configuration->layer_start_index[function_index];
        const int x_index = xx_start + i - layer_start_index;
        for(int j = 0; j < number_of_variables; j++){
            // Fill in ff
            ff[xx_start + i][j] = (Western_flux[i][j] - Eastern_flux[i][j])/(x_face[x_index+1]-x_face[x_index]) + Source[i][j];
        }
        // std::cout << Western_flux[i][0] << " " << Western_flux[i][2] << std::endl;
    }

    PetscCall(DMDAVecRestoreArrayDOF(da,local_x,&xx));
    PetscCall(DMDAVecRestoreArrayDOF(da,f,&ff));
    PetscCall(DMRestoreLocalVector(da, &local_x));

    // //////////////////////////////////////////////////////
    // Storing intermediate values so they are accesible for the output writer
    PetscScalar **fluxes, **sources;
    PetscCall(DMDAVecGetArrayDOF(da,Fluxes,&fluxes));
    PetscCall(DMDAVecGetArrayDOF(da,Source_term,&sources));
    for(int i = 0; i < xx_length; i++){
        for(int j = 0; j < number_of_variables; j++){
            fluxes[xx_start + i][j] = Western_flux[i][j];
            sources[xx_start + i][j] = Source[i][j];
        }
    }
    //?? not working yet
    // if(this->rank == comm_size-1){
    //     for(int j = 0; j < number_of_variables; j++){
    //         fluxes[xx_start+xx_length][j] = Eastern_flux[-1][j];
    //     }
    // }
    PetscCall(DMDAVecRestoreArrayDOF(da,Fluxes,&fluxes));
    PetscCall(DMDAVecRestoreArrayDOF(da,Source_term,&sources));
    
    
    PetscFunctionReturn(0);
}


PetscErrorCode RootFunction::calculate_jacobian(Vec x, Mat J, DM da) const{   
    const double eps = 1e-7; //?? move
    const double tiny = 1e-14; //?? move

    PetscInt number_of_variables, xx_length, xx_start, xx_global_length;

    PetscCall(DMDAGetInfo(da,PETSC_IGNORE,&xx_global_length,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,&number_of_variables,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE));
    PetscCall(DMDAGetCorners(da,&xx_start, NULL, NULL, &xx_length, NULL, NULL));
    // //Calculate the function vector for the original solution vector
    PetscCall(this->calculate_root_function(x,this->f, da)); //?? reuse function?
    int k = 0;

    for(int l = 0; l < comm_size; l++){
        for(int i = 0; i < nodes_per_process[l]; i++){
            for(int j = 0; j < number_of_variables; j++){
                PetscScalar xx_old, **xx, delta;
                PetscCall(DMDAVecGetArrayDOF(da,x,&xx));
                if(l == this->rank){
                    xx_old = xx[xx_start + i][j];
                    xx[xx_start+i][j] = xx[xx_start+i][j]*(1.0 + eps) + tiny;
                    delta = xx[xx_start+i][j] - xx_old;
                }
                PetscCall(DMDAVecRestoreArrayDOF(da,x,&xx));
                PetscCallMPI(MPI_Bcast(&delta,1,MPI_DOUBLE,l,PETSC_COMM_WORLD)); //Send delta to all other processes

                //Calculate deviated root function
                PetscCall(this->calculate_root_function(x,this->df, da)); 
                
                PetscCall(VecAXPBY(this->df,-1.0/delta,1.0/delta,this->f));

                const PetscScalar** dff;
                //Insert column k in the Jacobian
                PetscCall(DMDAVecGetArrayDOFRead(da,this->df,&dff));

                int idx; //Global row indices
                for(int m = xx_start; m < xx_start + xx_length; m++){
                    for(int n = 0; n < number_of_variables; n++){
                        idx = m*number_of_variables + n;
                        PetscCall(MatSetValue(J,idx,k,dff[m][n],INSERT_VALUES));
                    }
                }
                PetscCall(MatAssemblyBegin(J,MAT_FLUSH_ASSEMBLY));
                PetscCall(DMDAVecRestoreArrayDOFRead(da,this->df,&dff)); //Restore the array while waiting for the assembly
                PetscCall(MatAssemblyEnd(J,MAT_FLUSH_ASSEMBLY));

                //Restore solution vector
                PetscCall(DMDAVecGetArrayDOF(da,x,&xx));
                if(l == this->rank){
                    xx[xx_start + i][j] = xx_old; //Restore previous value
                }
                PetscCall(DMDAVecRestoreArrayDOF(da,x,&xx));
                k++;
            }

        }
    }

    PetscCall(MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY)); //Final assembly: compress matrix
    PetscCall(MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(0);
};

PetscErrorCode RootFunction::writer(Vec x, PetscScalar* xx_sol, PetscScalar* fluxes, PetscScalar* source, DM da) const{

    PetscInt number_of_variables, xx_length, xx_start, xx_global_length;
    PetscCall(DMDAGetInfo(da,PETSC_IGNORE,&xx_global_length,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,&number_of_variables,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE));
    PetscCall(DMDAGetCorners(da,&xx_start, NULL, NULL, &xx_length, NULL, NULL));
    const int electron_index = this->cell_configuration->electron_index;
    const int ion_index = this->cell_configuration->ion_index;

    if(this->rank > 0){
        PetscScalar **sol_vec;
        PetscCall(DMDAVecGetArrayDOFRead(da,x,&sol_vec));

        // // convert sol_Vec from 2D vector to 1D and send:
        PetscScalar xx_temp[nodes_per_process[this->rank]*number_of_variables];
        int k = 0;
        for(int i = 0; i < xx_length; i++){
            for(int j = 0; j < number_of_variables; j++){
                xx_temp[k] = sol_vec[xx_start+i][j];
                k++;
            }
        }

        PetscCallMPI(MPI_Send(xx_temp,nodes_per_process[this->rank]*number_of_variables,MPI_DOUBLE,0,0,PETSC_COMM_WORLD));
        
        PetscCall(DMDAVecRestoreArrayDOFRead(da,x,&sol_vec));
    }

    if(this->rank == 0){

        PetscScalar **sol_vec;
        PetscCall(DMDAVecGetArrayDOFRead(da,x,&sol_vec));

        for(int i = 0; i < xx_length; i++){
            for(int j = 0; j < number_of_variables; j++){
                xx_sol[i*number_of_variables + j] = sol_vec[i][j];
            }
        }

        int starting_node = nodes_per_process[0];
        for(int comm = 1; comm < comm_size; comm++){
            int k = 0;
            PetscScalar xx_temp[nodes_per_process[comm]*number_of_variables];
            // MPI recieve commands and store in 1D vector
            PetscCallMPI(MPI_Recv(&xx_temp[0],nodes_per_process[comm]*number_of_variables,MPI_DOUBLE,comm,0,PETSC_COMM_WORLD,MPI_STATUS_IGNORE));
            for(int i = starting_node; i < starting_node + nodes_per_process[comm]; i++){
                for(int j = 0; j < number_of_variables; j++){
                    xx_sol[i*number_of_variables + j] = xx_temp[k];
                    k++;
                }
            }
            starting_node += nodes_per_process[comm];
        }
        //Iterate over all nodes
        for(int i = 0; i < xx_global_length; i++){
            const int function_index = node_layer_index[i];
            const int layer_start_index = this->cell_configuration->layer_start_index[function_index];
            //Calculate western fluxes
            for(int j = 0; j < number_of_variables; j++){
                fluxes[i*number_of_variables + j] = 0.0;
                source[i*number_of_variables + j] = 0.0; // Set all values to zero
            }
            // //The flux is calculated in the positive x direction (J = - D dc/dx)
            if(i == 0){
                this->root_sub_functions[function_index]->calculate_flux(&xx_sol[i*number_of_variables],i - layer_start_index, &fluxes[i*number_of_variables],electron_index,ion_index);
            }
            else{
                this->root_sub_functions[function_index]->calculate_flux(   &xx_sol[i*number_of_variables],&xx_sol[(i-1)*number_of_variables], i - layer_start_index,
                                                                            &fluxes[i*number_of_variables],electron_index,ion_index);
            }
            this->root_sub_functions[function_index]->calculate_source_term(&xx_sol[i*number_of_variables], &source[i*number_of_variables], electron_index, ion_index);
            std::cout << fluxes[i*number_of_variables] << " " << fluxes[i*number_of_variables+1] << " " << fluxes[i*number_of_variables+2] << " " << fluxes[i*number_of_variables+3] << " " << std::endl;
        }

    }
    PetscFunctionReturn(0);
}

const std::string RootFunction::get_petsc_convergence_reasons(const int key) const{
    std::string output;
    if(this->petsc_convergence_reasons.count(key) == 1){ //Check if key is present in petsc_convergence_reasons
        output = this->petsc_convergence_reasons.at(key);
    }
    else{
        std::cerr << "Error: value requested for " << key << " not present in unordered map (get_petsc_convergence_reasons)" << std::endl;
    }

    return output;
};

typedef struct{
    DM da;
    std::unique_ptr<RootFunction> rf_wrap;
}   AppCtx;


//?? make wrappers more efficient?

//Wrapper around calculate_root_function and calculate_jacobian to ensure compatibility with PetSc for class member functions
//NOT a member of the RootFunction class itself
PetscErrorCode root_function_wrapper(SNES snes_context, Vec x, Vec f, void * obj){

    AppCtx     *Ctx = (AppCtx *)obj;
    PetscCall((Ctx->rf_wrap)->calculate_root_function(x,f,Ctx->da));

    PetscFunctionReturn(0);
}

PetscErrorCode jacobian_wrapper(SNES snes_context, Vec x, Mat J, Mat PC, void * obj){

    AppCtx     *Ctx = (AppCtx *)obj;
    PetscCall((Ctx->rf_wrap)->calculate_jacobian(x,J,Ctx->da));

    PetscFunctionReturn(0);
};

PetscErrorCode write_wrapper(const int rank, const std::string& input_folder, const std::string& output_folder, Vec x, void * obj){

    AppCtx     *Ctx = (AppCtx *)obj;
    if(rank == 0){
        //Check for existence of output directory
        if(std::filesystem::is_directory(output_folder)){
            std::filesystem::remove_all(output_folder);
            std::filesystem::create_directory(output_folder);
        }
        else{
            std::filesystem::create_directory(output_folder);
        }

        //Copy input folder into output
        std::filesystem::copy(input_folder,output_folder + "input/");

    }
    PetscInt number_of_variables, xx_length, xx_start, xx_global_length;
    PetscCall(DMDAGetInfo(Ctx->da,PETSC_IGNORE,&xx_global_length,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,&number_of_variables,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE));
    PetscScalar xx[xx_global_length*number_of_variables];
    PetscScalar ff[xx_global_length*number_of_variables];
    PetscScalar ss[xx_global_length*number_of_variables];
    PetscCall((Ctx->rf_wrap)->writer(x,xx,ff,ss,Ctx->da));

    PetscFunctionReturn(0);
}
