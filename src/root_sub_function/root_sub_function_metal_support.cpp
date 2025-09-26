#include "root_sub_function_metal_support.h"

RootSubFunctionMetalSupport::RootSubFunctionMetalSupport( const int _rank, const std::vector<double>& _x_center, const std::vector<double>& _x_face,
                                                          const std::vector<std::string>& gas_species, const double temperature, const bool _galvanostatic,
                                                          const double _applied_I_or_V, const double area,
                                                          const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration):
RootSubFunctionPorous(gas_species,temperature,area,layer_index,cell_configuration), 
rank(_rank),x_center(_x_center),x_face(_x_face),galvanostatic(_galvanostatic),applied_I_or_V(_applied_I_or_V){

    //Derive fields from constructor input
    this->nodes = x_center.size();
    this->number_of_gas_species = this->gas_species.size();

    //Determine whether the preceding layer is dense
    if(layer_index > 0){
        this->preceding_layer_is_dense = cell_configuration->get_layer_is_dense(layer_index-1);
    }

    //Determine whether the preceding layer conducts electrons. If so, determine the index
    if(layer_index > 0){
        const int preceding_layer_number_of_gas_species = cell_configuration->get_layer_number_of_gas_species(layer_index-1);
        const std::pair<int,int> preceding_layer_charge_carriers = cell_configuration->get_layer_charge_carriers(layer_index-1);
        if(preceding_layer_charge_carriers.first > 0){
            this->preceding_layer_conducts_electrons = 1;
            this->preceding_layer_electron_index = preceding_layer_number_of_gas_species;
            //Note that if there are multiple electron conductors in the preceding layer, the first is taken (which is considered to 
            //correspond to the main electron conductor)            
        }
        else{this->preceding_layer_conducts_electrons = 0;}
    }

    this->set_effective_binary_diffusion_coefficients_constant_quotient(this->porosity,this->tortuosity);

};

void RootSubFunctionMetalSupport::calculate_flux(   const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, 
                                                    PetscScalar* flux, const int electron_index, const int ion_index){
    

    if(x_index > 0){ //Fluxes throughout the layer
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
        const double sigma_e_w = (sigma_e_P + sigma_e_W)/2.0; //Take average
        flux[electron_index] = -(sigma_e_w/(this->z[0]*this->F))*
                                            ((mu_e_P - mu_e_W)/(this->x_center[x_index] - this->x_center[x_index-1]));

        flux[ion_index] = center_node[ion_index]-west_node[ion_index];
    }
    else{ //Flux at the interface between the porous support and the previous layer
        if(this->preceding_layer_is_dense){
            std::cerr << "Attempt to calculate fluxes between dense layer and porous metal support. Insert electro-active layer in between." << std::endl;
        }
        else{//Porous-porous interface
            const PetscScalar* c_i_P = center_node;
            PetscScalar c_i_w[this->number_of_gas_species], nn[this->number_of_gas_species];
            const PetscScalar* c_i_W = west_node;
            const double Delta_x_wW = this->preceding_layer->get_final_face_to_center_width();
            
            this->calculate_concentration_face(c_i_W,c_i_P,c_i_w);
            this->calculate_gas_flux(c_i_P,c_i_W,c_i_w,this->x_center[x_index],this->x_face[x_index]-Delta_x_wW, nn);
            for(int i = 0; i < this->number_of_gas_species; i++){
                flux[i] = nn[i];
            }         
            
            //This implementation is the same as for the case of oxygen anions in YSZ and the preceding layer.
            //Electrons
            const PetscScalar mu_e_W = west_node[electron_index];
            const PetscScalar mu_e_P = center_node[electron_index];
            const double sigma_e_P = this->calculate_electron_conductivity(center_node);
            const double sigma_e_W = this->preceding_layer->calculate_electron_conductivity(west_node);
            const double Delta_x_Pw = this->x_center[0] - this->x_face[0];
            const double quotient = (sigma_e_P*Delta_x_wW)/(sigma_e_W*Delta_x_Pw);

            flux[electron_index] =  -((sigma_e_P)/(this->z[0]*this->F))*
                                                (mu_e_P - (quotient*mu_e_P + mu_e_W)/(1.0+quotient))/(Delta_x_Pw);
            
            flux[ion_index] = center_node[ion_index]-0.0;
        }
    }
};

void RootSubFunctionMetalSupport::calculate_flux(   const PetscScalar* center_node, const int x_index, PetscScalar* flux,
                                                    const int electron_index, const int ion_index){
    //This function overload is used for the last node (i.e. rank = comm size - 1, node_index = final node and x_index = nodes - 1
    //- the interface with the fuel compartment) to calculate the eastern/final flux 
    
    //Note that at the fuel compartment interface the electrochemical electron potential is always set to zero. At the air compartement,
    //either the electrochemical potential is fixed to E_cell (potentiostatic mode), or the electron flux is set to the applied current
    //(galvanostatic mode) 

    //Extract values
    const PetscScalar mu_e_P = center_node[electron_index];

    if(x_index == this->nodes - 1){ //Interface with fuel compartment
        //Mass flux
        const double inlet_concentration = this->gas_compartment->inlet_pressure/(this->R*this->temperature);
        double c_i_bc[this->number_of_gas_species];
        for(int i = 0; i < this->number_of_gas_species; i++){
            c_i_bc[i] = this->gas_compartment->inlet_mole_fractions[i]*inlet_concentration;
        }
        const PetscScalar* c_i_P = center_node;
        //Mass flux - solve the DGM linear system //?? copied from LSCF
        PetscScalar c_i_e[this->number_of_gas_species], nn[this->number_of_gas_species];
        this->calculate_concentration_face(c_i_P,c_i_bc,c_i_e);
        this->calculate_gas_flux(c_i_bc,c_i_P,c_i_e,this->x_face[x_index+1],this->x_center[x_index], nn);
        for(int i = 0; i < this->number_of_gas_species; i++){
            flux[i] = nn[i];
        }

        

        //Electron flux (equal for potentiostatic and galvanostatic)
        const double sigma_e_P = this->calculate_electron_conductivity(center_node); //?? sigma_e_e can technically be calculated (interface values from gas compartment)
        flux[electron_index] = -(sigma_e_P/(this->z[0]*this->F))*
                                            ((0.0 - mu_e_P)/(this->x_face[x_index+1] - this->x_center[x_index]));

    }
    else{
        std::cerr << "Warning: calculate_flux function overload for the gas compartment interfaces called for invalid x_index" << std::endl;
    }     
}

void RootSubFunctionMetalSupport::calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const{
    
}

const double RootSubFunctionMetalSupport::get_final_face_to_center_width() const{
    return this->x_face[this->nodes] - this->x_center[this->nodes-1];
};

const double RootSubFunctionMetalSupport::calculate_electron_conductivity(const PetscScalar* node) const{
    double sigma0 = 3.27e6 - 1065.3*this->temperature;
    double sigma = 1e2*sigma0*pow((1 - this->porosity), 1.5);
    
    //?? implement no species situation
    return sigma;
};