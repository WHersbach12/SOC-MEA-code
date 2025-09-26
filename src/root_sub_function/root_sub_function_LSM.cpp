#include "root_sub_function_LSM.h"

RootSubFunctionLSM::RootSubFunctionLSM(   const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
                                            const std::vector<std::string>& gas_species, const double temperature, const double area, const bool _galvanostatic,
                                            const double _applied_I_or_V, const int layer_index, 
                                            const std::shared_ptr<CellConfiguration>& cell_configuration):
RootSubFunctionPorous(gas_species,temperature,area,layer_index,cell_configuration),
rank(_rank),x_center(_x_center),x_face(_x_face),galvanostatic(_galvanostatic),applied_I_or_V(_applied_I_or_V){
    
    //Derive fields from constructor input
    this->nodes = x_center.size();


    if(layer_index > 0){
        this->preceding_layer_is_dense = cell_configuration->get_layer_is_dense(layer_index-1);
    }

    this->set_effective_binary_diffusion_coefficients_constant_quotient(this->porosity,this->tortuosity);
};

RootSubFunctionLSM::~RootSubFunctionLSM(){};    

void RootSubFunctionLSM::calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, PetscScalar* flux, 
                    const int electron_index, const int ion_index){
    //Calculate the flux at the western cell face, an array with length equal to center_node_length
    
    //x_index: x index of cell center P

    if(x_index > 0){ //Flux at all cell faces inside LSCF
        //Extract values
        const PetscScalar* c_i_P = center_node;
        const PetscScalar mu_e_P = center_node[electron_index];
        const PetscScalar* c_i_W = west_node;
        const PetscScalar mu_e_W = west_node[electron_index];

        //Mass flux - solve the DGM linear system //?? copied from LSCF
        PetscScalar c_i_w[this->number_of_gas_species], nn[this->number_of_gas_species];
        this->calculate_concentration_face(c_i_W,c_i_P,c_i_w);
        this->calculate_gas_flux(c_i_P,c_i_W,c_i_w,this->x_center[x_index],this->x_center[x_index-1], nn);
        for(int i = 0; i < this->number_of_gas_species; i++){
            flux[i] = nn[i];
        }
    
        //Electron flux
        const double sigma_e_P = this->calculate_electron_conductivity(center_node);
        const double sigma_e_W = this->calculate_electron_conductivity(west_node);
        const double sigma_e_w = (sigma_e_P + sigma_e_W)/2.0; //Take average //?? take average beforehand?
        flux[electron_index] = -(sigma_e_w/(this->z[0]*this->F))*
                                             ((mu_e_P - mu_e_W)/(this->x_center[x_index] - this->x_center[x_index-1]));

        //Oxygen anion flux
        const PetscScalar mu_o_P = center_node[ion_index];
        const PetscScalar mu_o_W = west_node[ion_index];
        //Oxygen anion flux
        flux[ion_index] = mu_o_P - mu_o_W;

    }
    else{ //As this version of the overload is called, there must exist a preceding node and LSCF does not interface the air compartment
        if(this->preceding_layer_is_dense){  
            //This implementation is the same as for the case of oxygen anions in YSZ and the preceding layer.
            //Electrons
            const PetscScalar mu_e_W = west_node[electron_index];
            const PetscScalar mu_e_P = center_node[electron_index];
            const double sigma_e_P = this->calculate_electron_conductivity(center_node);
            const double sigma_e_W = this->preceding_layer->calculate_electron_conductivity(west_node);
            const double Delta_x_Pw = this->x_center[x_index] - this->x_face[x_index];
            const double Delta_x_wW = this->preceding_layer->get_final_face_to_center_width();
            const double quotient = (sigma_e_P*Delta_x_wW)/(sigma_e_W*Delta_x_Pw);

            flux[electron_index] =  -((sigma_e_P)/(this->z[0]*this->F))*
                                                (mu_e_P - (quotient*mu_e_P + mu_e_W)/(1.0+quotient))/(Delta_x_Pw);

            //Oxygen anions
            const PetscScalar mu_o_P = center_node[ion_index];
            const PetscScalar mu_o_W = west_node[ion_index];
            //Oxygen anion flux
            flux[ion_index] = mu_o_P - 0.0;    
        }
        else{
            const PetscScalar* c_i_P = center_node;
            PetscScalar c_i_w[this->number_of_gas_species], nn[this->number_of_gas_species];
            const PetscScalar* c_i_W = west_node;
            const double Delta_x_wW = this->preceding_layer->get_final_face_to_center_width();
            
            this->calculate_concentration_face(c_i_W,c_i_P,c_i_w);
            this->calculate_gas_flux(c_i_P,c_i_W,c_i_w,this->x_center[x_index],this->x_face[x_index]-Delta_x_wW, nn);
            // this->calculate_concentration_face(c_i_w,c_i_W,c_i_w_av);
            // this->preceding_layer->calculate_gas_flux(c_i_w,c_i_W,c_i_w_av,this->x_face[x_index],this->x_face[x_index]-Delta_x_wW, nn_2);
            for(int i = 0; i < this->number_of_gas_species; i++){
                flux[i] = nn[i];
                // flux[this->number_of_gas_species+i] = nn_2[i] - nn[i];s
                // std::cout << c_i_w[i]  - c_i_W[i]<< std::endl;
            }
            
            //This implementation is the same as for the case of oxygen anions in YSZ and the preceding layer.
            //Electrons
            const PetscScalar mu_e_W = west_node[electron_index];
            const PetscScalar mu_e_P = center_node[electron_index];
            const double sigma_e_P = this->calculate_electron_conductivity(center_node);
            const double sigma_e_W = this->preceding_layer->calculate_electron_conductivity(west_node);
            const double Delta_x_Pw = this->x_center[x_index] - this->x_face[x_index];
            const double quotient = (sigma_e_P*Delta_x_wW)/(sigma_e_W*Delta_x_Pw);

            flux[electron_index] =  -((sigma_e_P)/(this->z[0]*this->F))*
                                                (mu_e_P - (quotient*mu_e_P + mu_e_W)/(1.0+quotient))/(Delta_x_Pw);

            //Oxygen anions
            const PetscScalar mu_o_P = center_node[ion_index];
            const PetscScalar mu_o_W = west_node[ion_index];
            //Oxygen anion flux
            flux[ion_index] = mu_o_P - 0.0;
        }
    }


};

void RootSubFunctionLSM::calculate_flux(const PetscScalar* center_node, const int x_index, PetscScalar* flux, const int electron_index, const int ion_index){
    //This function overload is used when no preceding node exists (i.e. rank = 0, node_index = 0 and x_index = 0
    //- the interface with the air compartment) to calculate the western flux or when no preceding node is required 
    //to calculate the eastern (!) flux (i.e. rank = comm_size - 1, node_index = final node and x_index = nodes - 1 
    //- the interface with the fuel compartment)
    
    //Note that at the fuel compartment interface the electrochemical electron potential is always set to zero. At the air compartement,
    //either the electrochemical potential is fixed to E_cell (potentiostatic mode), or the electron flux is set to the applied current
    //(galvanostatic mode) 

    //Extract values
    const PetscScalar mu_e_P = center_node[electron_index];
    const PetscScalar mu_o_P = center_node[ion_index];

    if(x_index == 0){ //Interface with air compartment
        //Mass flux
        const double inlet_concentration = this->gas_compartment->inlet_pressure/(this->R*this->temperature);
        PetscScalar c_i_bc[this->number_of_gas_species];
        for(int i = 0; i < this->number_of_gas_species; i++){
            c_i_bc[i] = (this->gas_compartment->inlet_mole_fractions[i])*inlet_concentration;
        }
        const PetscScalar* c_i_P = center_node;
        //Mass flux - solve the DGM linear system //?? copied from LSCF
        PetscScalar c_i_w[this->number_of_gas_species], nn[this->number_of_gas_species];
        this->calculate_concentration_face(c_i_bc,c_i_P,c_i_w);
        this->calculate_gas_flux(c_i_P,c_i_bc,c_i_w,this->x_center[x_index],this->x_face[x_index], nn);
        for(int i = 0; i < this->number_of_gas_species; i++){
            flux[i] = nn[i];
        }
    

        //Electron flux
        if(this->galvanostatic){ //Applied current
            flux[electron_index] = this->applied_I_or_V;
        }
        else{ //Applied potential
            const double sigma_e_P = this->calculate_electron_conductivity(center_node); //?? sigma_e_w can technically be calculated
            const double E_cell = this->applied_I_or_V*this->F;

            flux[electron_index] = -(sigma_e_P/(this->z[0]*this->F))*
                                                ((mu_e_P - E_cell)/(this->x_center[x_index] - this->x_face[x_index]));
        }

        //Oxygen anion flux
        flux[ion_index] = mu_o_P - 0.0;

    }
    else if(x_index == this->nodes - 1){ //Interface with fuel compartment
        //Mass flux

        //Electron flux (equal for potentiostatic and galvanostatic)
        const double sigma_e_P = this->calculate_electron_conductivity(center_node); //?? sigma_e_e can technically be calculated (interface value from gas compartment)
        flux[electron_index] = -(sigma_e_P/(this->z[0]*this->F))*
                                            ((0.0 - mu_e_P)/(this->x_face[x_index+1] - this->x_center[x_index]));

        //Oxygen anion flux
        flux[ion_index] = 0.0;
    }
    else{
        std::cerr << "Warning: calculate_flux function overload for the gas compartment interfaces called for invalid x_index" << std::endl;
    }
}; 

const double RootSubFunctionLSM::get_final_face_to_center_width() const{
    return this->x_face[this->nodes] - this->x_center[this->nodes-1];
};

void RootSubFunctionLSM::calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const{


};


const double RootSubFunctionLSM::calculate_electron_conductivity(const PetscScalar* node) const{

    // [Zhu 2008] Modelling distributed charge-transfer proceses in SOFC membrane electrode assemblies
    double sigma0 = 8.855e5/this->temperature * exp(-1082.5/this->temperature);
    double sigma = 1e2*sigma0*pow((1 - this->porosity), 1.5);

    return sigma;
};
