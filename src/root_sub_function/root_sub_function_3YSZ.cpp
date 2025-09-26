#include "root_sub_function_3YSZ.h"

RootSubFunction3YSZ::RootSubFunction3YSZ( const int _rank, const std::vector<double>& _x_center, const std::vector<double>& _x_face, const double _temperature,
                                        const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration):
rank(_rank),x_center(_x_center),x_face(_x_face),temperature(_temperature){

    //Derive fields from constructor input
    this->nodes = x_center.size();

    //Identify the index for the oxygen anion electrochemical potential in the unit of the preceding layer
    if(layer_index > 0){
        const int preceding_layer_number_of_gas_species = cell_configuration->get_layer_number_of_gas_species(layer_index-1);
        const std::pair<int,int> preceding_layer_charge_carriers = cell_configuration->get_layer_charge_carriers(layer_index-1);
        this->preceding_layer_oxygen_anion_index = preceding_layer_number_of_gas_species + preceding_layer_charge_carriers.first;
        //Note that if there are multiple oxygen anion conductors in the preceding layer, the first is taken (which is considered to 
        //correspond to the main oxygen anion conductor)
    }
    

};


void RootSubFunction3YSZ::calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node,
                                        const int x_index, PetscScalar* flux, const int electron_index, const int ion_index){
    
    if(x_index > 0){ //Flux at all cell faces inside YSZ
        //Extract values
        const PetscScalar mu_o_P = center_node[ion_index];
        const PetscScalar mu_o_W = west_node[ion_index];

        //Oxygen anion flux
        const double sigma_o_P = this->calculate_ion_conductivity(center_node);
        const double sigma_o_W = this->calculate_ion_conductivity(west_node);
        const double sigma_o_w = (sigma_o_P + sigma_o_W)/2.0; //Take average
        flux[ion_index] =   -(sigma_o_w/(this->z[0]*this->F))*
                    ((mu_o_P - mu_o_W)/(this->x_center[x_index] - this->x_center[x_index-1]));

        for(int i =0; i<ion_index; i++){
            flux[i] = center_node[i]-west_node[i];
        }

    }
    else{ //Flux at the interface between YSZ and the preceding layer
        // other|YSZ
        //  | * | * |
        //    W w P e
        //      ^
        //      Interface
   
        //Extract values
        const PetscScalar mu_o_P = center_node[ion_index];
        const PetscScalar mu_o_W = west_node[ion_index];

        //Calculate conductivity
        const double sigma_o_P = this->calculate_ion_conductivity(center_node);
        //It is assumed that sigma_o_P equals sigma_o_w, because the conductivity is assumed to
        //change little from P to w
        const double sigma_o_W = this->preceding_layer->calculate_ion_conductivity(west_node);
        //It is assumed that sigma_o_W equals sigma_o_w, which is valid as long as sigma_o_W
        //is not a function of mu_o_W (YSZ is a dense oxygen anion conductor, a Neumann boundary
        //condition applies (slope = 0) for all other dependent variables in the preceding layer
        //and thus W ~ w)

        //Calculate Delta_x
        const double Delta_x_Pw = this->x_center[x_index] - this->x_face[x_index];
        const double Delta_x_wW = this->preceding_layer->get_final_face_to_center_width();
        const double quotient = (sigma_o_P*Delta_x_wW)/(sigma_o_W*Delta_x_Pw);

        flux[ion_index] =   -((sigma_o_P)/(this->z[0]*this->F))*
                    (mu_o_P - (quotient*mu_o_P + mu_o_W)/(1.0+quotient))/(Delta_x_Pw);

        for(int i =0; i<ion_index; i++){
            flux[i] = center_node[i]-0.0;
        }
    }

};


void RootSubFunction3YSZ::calculate_flux(const PetscScalar* center_node, const int x_index, PetscScalar* flux, const int electron_index, const int ion_index){

    //?? Not configured yet
    std::cerr << "Attempt to calculate the flux at the western interface between YSZ and the air or fuel compartment - not yet configured" << std::endl;

}; 

const double RootSubFunction3YSZ::get_final_face_to_center_width() const{
    return this->x_face[this->nodes] - this->x_center[this->nodes-1];
};


void RootSubFunction3YSZ::calculate_source_term(const PetscScalar* unit, PetscScalar* source, const int electron_index, const int ion_index) const{
    //No source term in YSZ

};

const double RootSubFunction3YSZ::calculate_ion_conductivity(const PetscScalar* unit) const{

    // [Lay-Gindler 2013] https://doi.org/10.1016/j.ijhydene.2013.03.162
    double sigma = 45e2 * exp(-69.837e3/this->R/this->temperature);
       
    return sigma;
};


