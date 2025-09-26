#include "output_writer.h"

OutputWriter::OutputWriter(const int _rank, const std::string& _output_folder, const std::string& input_folder):
rank(_rank),output_folder(_output_folder){


    if(this->rank == 0){
        //Check for existence of output directory
        if(std::filesystem::is_directory(this->output_folder)){
            std::filesystem::remove_all(this->output_folder);
            std::filesystem::create_directory(this->output_folder);
        }
        else{
            std::filesystem::create_directory(this->output_folder);
        }

        //Copy input folder into output
        std::filesystem::copy(input_folder,output_folder + "input/");

    }
};


PetscErrorCode OutputWriter::write_solution_to_csv( const int* nodes_per_process, const int comm_size, const int solution_vector_length, 
                                                    const std::shared_ptr<CellConfiguration>& cell_configuration, Vec x, Vec Fluxes, Vec Source_term, DM da) const{

    PetscInt number_of_variables, xx_length, xx_start, xx_global_length;
    PetscCall(DMDAGetInfo(da,PETSC_IGNORE,&xx_global_length,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,&number_of_variables,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE));
    PetscCall(DMDAGetCorners(da,&xx_start, NULL, NULL, &xx_length, NULL, NULL));

    if(this->rank > 0){
        PetscScalar **sol_vec, **fluxes, **sources;
        PetscCall(DMDAVecGetArrayDOFRead(da,x,&sol_vec));
        PetscCall(DMDAVecGetArrayDOFRead(da, Fluxes, &fluxes));
        PetscCall(DMDAVecGetArrayDOFRead(da, Source_term, &sources));

        // // convert sol_Vec from 2D vector to 1D and send:
        PetscScalar xx_temp[nodes_per_process[this->rank]*number_of_variables];
        PetscScalar ff_temp[nodes_per_process[this->rank]*number_of_variables];
        PetscScalar ss_temp[nodes_per_process[this->rank]*number_of_variables];
        int k = 0;
        for(int i = 0; i < xx_length; i++){
            for(int j = 0; j < number_of_variables; j++){
                xx_temp[k] = sol_vec[xx_start+i][j];
                ff_temp[k] = fluxes[xx_start+i][j];
                ss_temp[k] = sources[xx_start+i][j];
                k++;
            }
        }

        PetscCallMPI(MPI_Send(xx_temp,nodes_per_process[this->rank]*number_of_variables,MPI_DOUBLE,0,0,PETSC_COMM_WORLD));
        PetscCallMPI(MPI_Send(ff_temp,nodes_per_process[this->rank]*number_of_variables,MPI_DOUBLE,0,0,PETSC_COMM_WORLD));
        PetscCallMPI(MPI_Send(ss_temp,nodes_per_process[this->rank]*number_of_variables,MPI_DOUBLE,0,0,PETSC_COMM_WORLD));
        
        PetscCall(DMDAVecRestoreArrayDOFRead(da,x,&sol_vec));
        PetscCall(DMDAVecRestoreArrayDOFRead(da,Fluxes,&fluxes));
        PetscCall(DMDAVecRestoreArrayDOFRead(da,Source_term,&sources));
    }

    if(this->rank == 0){
        PetscScalar **sol_vec, **fluxes, **sources;
        PetscScalar xx[xx_global_length][number_of_variables], ff[xx_global_length][number_of_variables], ss[xx_global_length][number_of_variables];
        PetscCall(DMDAVecGetArrayDOFRead(da,x,&sol_vec));
        PetscCall(DMDAVecGetArrayDOFRead(da, Fluxes, &fluxes));
        PetscCall(DMDAVecGetArrayDOFRead(da, Source_term, &sources));
        for(int i = 0; i < xx_length; i++){
            for(int j = 0; j < number_of_variables; j++){
                xx[i][j] = sol_vec[i][j];
                ff[i][j] = fluxes[i][j];
                ss[i][j] = sources[i][j];
            }
        }
        int starting_node = nodes_per_process[0];
        for(int comm = 1; comm < comm_size; comm++){
            int k = 0;
            PetscScalar xx_temp[nodes_per_process[comm]*number_of_variables];
            PetscScalar ff_temp[nodes_per_process[comm]*number_of_variables];
            PetscScalar ss_temp[nodes_per_process[comm]*number_of_variables];
            // MPI recieve commands and store in 1D vector
            PetscCallMPI(MPI_Recv(&xx_temp[0],nodes_per_process[comm]*number_of_variables,MPI_DOUBLE,comm,0,PETSC_COMM_WORLD,MPI_STATUS_IGNORE));
            PetscCallMPI(MPI_Recv(&ff_temp[0],nodes_per_process[comm]*number_of_variables,MPI_DOUBLE,comm,0,PETSC_COMM_WORLD,MPI_STATUS_IGNORE));
            PetscCallMPI(MPI_Recv(&ss_temp[0],nodes_per_process[comm]*number_of_variables,MPI_DOUBLE,comm,0,PETSC_COMM_WORLD,MPI_STATUS_IGNORE));
            for(int i = starting_node; i < starting_node + nodes_per_process[comm]; i++){
                for(int j = 0; j < number_of_variables; j++){
                    xx[i][j] = xx_temp[k];
                    ff[i][j] = ff_temp[k];
                    ss[i][j] = ss_temp[k];
                    k++;
                }
            }
            starting_node += nodes_per_process[comm];
        }
        //Create array for all central x-positions, to be owned by rank 0
        double xc[solution_vector_length];
        double xf[solution_vector_length+1];
        std::cout << "Cell voltage is " << xx[0][2]/96485.33 << " V at current density of " << ff[0][2]/1e4 << std::endl;
        const int layer_count = cell_configuration->get_layer_count();
        for(int i = 0; i < layer_count; i++){ //Iterate over all layers
            const int layer_nodes = cell_configuration->get_layer_nodes(i);
            const std::vector<double> x_center = cell_configuration->get_layer_x_center(i);
            const std::vector<double> x_face = cell_configuration->get_layer_x_face(i);
            const int layer_start_index = cell_configuration->layer_start_index[i];
            for(int j = 0; j < layer_nodes; j++){ //Iterate over all nodes
                    xc[layer_start_index+ j] = x_center[j];
                    xf[layer_start_index+ j] = x_face[j];
            }
            xf[-1] = x_face[-1];
        }
        const std::string file_name = "solution_vector.txt";
        std::ofstream out(this->output_folder + file_name); //Open file
        out << std::setprecision(this->precision); //Set precision for writing values
        if (out.is_open()){
            for(int i = 0; i < solution_vector_length; i++){ //Iterate over all lines
                out << xc[i];
                for(int j = 0; j < number_of_variables; j++){
                    out << ' '<< xx[i][j] ;
                }
                for(int j = 0; j < number_of_variables; j++){
                    out << ' '<< ss[i][j] ;
                }
                out << ' ';
                out << xf[i];
                for(int j = 0; j < number_of_variables; j++){
                    out << ' '<< ff[i][j] ;
                }   
                out << std::endl;
            }
            out.close();
        }
        else{
            std::cerr << "Warning: unable to open output file (" << this->output_folder << file_name << ")" << std::endl;
        }
        PetscCall(DMDAVecRestoreArrayDOFRead(da,x,&sol_vec));
        PetscCall(DMDAVecRestoreArrayDOFRead(da,Fluxes,&fluxes));
        PetscCall(DMDAVecRestoreArrayDOFRead(da,Source_term,&sources));
    }

    PetscFunctionReturn(0);
};