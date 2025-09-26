#include "root_sub_function_CGO.h"

RootSubFunctionCGO::RootSubFunctionCGO(   const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
                                            const std::vector<std::string>& gas_species, const double temperature, const double area, const bool _galvanostatic,
                                            const double _applied_I_or_V, const int layer_index, 
                                            const std::shared_ptr<CellConfiguration>& cell_configuration):
RootSubFunctionPorous(gas_species,temperature,area,layer_index,cell_configuration),
rank(_rank),x_center(_x_center),x_face(_x_face),galvanostatic(_galvanostatic),applied_I_or_V(_applied_I_or_V){
    
    //Derive fields from constructor input
    this->nodes = x_center.size();

    //Identify index for oxygen species
    for(int i = 0; i < this->number_of_gas_species; i++){
        if(this->gas_species[i] == "O2"){
            this->oxygen_index = i;
            break;
        }
    }
    if(this->oxygen_index < 0){
        std::cerr << "Warning: oxygen species not found for CGO - not compatible with electron conductivity correlation (RootSubFunctionCGO)" << std::endl;
    }

    //Determine whether the preceding layer is dense
    if(layer_index > 0){
        this->preceding_layer_is_dense = cell_configuration->get_layer_is_dense(layer_index-1);
    }

    //Determine whether the preceding layer conducts electrons, oxygen anions. If so, determine the index
    const int preceding_layer_number_of_gas_species = cell_configuration->get_layer_number_of_gas_species(layer_index-1);
    if(layer_index > 0){
        const std::pair<int,int> preceding_layer_charge_carriers = cell_configuration->get_layer_charge_carriers(layer_index-1);
        if(preceding_layer_charge_carriers.first > 0){
            this->preceding_layer_conducts_electrons = 1;
            this->preceding_layer_electron_index = preceding_layer_number_of_gas_species;
            //Note that if there are multiple electron conductors in the preceding layer, the first is taken (which is considered to 
            //correspond to the main electron conductor)            
        }
        else{this->preceding_layer_conducts_electrons = 0;}
        if(preceding_layer_charge_carriers.second > 0){
            this->preceding_layer_conducts_oxygen_anions = 1;
            this->preceding_layer_oxygen_anion_index = preceding_layer_number_of_gas_species + preceding_layer_charge_carriers.first;
            //Note that if there are multiple oxygen anion conductors in the preceding layer, the first is taken (which is considered to 
            //correspond to the main oxygen anion conductor)
        }
        else{this->preceding_layer_conducts_oxygen_anions = 0;}
    }

    //Percolation probability based on [Bertei 2011] A comparative study and an extended theory of percolation for random packings of rigid spheres
    // and [Bertei 2011] Percolation theory in SOFC composite electrodes
    const double nr_particles = (1-this->porosity)/(3.1416/3*this->particle_diameter*this->particle_diameter);
    const double particle_area = 6*(1-this->porosity)/this->particle_diameter;//3.1416*nr_particles*this->particle_diameter*this->particle_diameter;

};

RootSubFunctionCGO::~RootSubFunctionCGO(){};    

void RootSubFunctionCGO::calculate_flux(   const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, 
                                            PetscScalar* flux, const int electron_index, const int ion_index){
    //Calculate the flux at the western cell face, an array with length equal to center_node_length
    
    //x_index: x index of cell center P

    if(x_index > 0){ //Flux at all cell faces inside CGO
        //Extract values
        const PetscScalar* c_i_P = center_node;
        const PetscScalar mu_e_P = center_node[electron_index];
        const PetscScalar mu_o_P = center_node[ion_index];
        const PetscScalar* c_i_W = west_node;
        const PetscScalar mu_e_W = west_node[electron_index];
        const PetscScalar mu_o_W = west_node[ion_index];

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
        const double sigma_o_P = this->calculate_ion_conductivity(center_node);
        const double sigma_o_W = this->calculate_ion_conductivity(west_node);
        const double sigma_o_w = (sigma_o_P + sigma_o_W)/2.0; //Take average
        flux[ion_index] = -(sigma_o_w/(this->z[1]*this->F))*
                                             ((mu_o_P - mu_o_W)/(this->x_center[x_index] - this->x_center[x_index-1]));

    }
    else{ //As this version of the overload is called, there must exist a preceding node and CGO does not interface the air compartment
        if(this->preceding_layer_is_dense){ //Dense-porous interface
            //Mass
            for(int i = 0; i < this->number_of_gas_species; i++){
                flux[i] = 0.0; //No gas-species flux through the interface
            }

            //Electrons
            if(this->preceding_layer_conducts_electrons){
                std::cerr << "Attempt to calculate electron the flux at the western interface between an electron conducting phase and Ni-YSZ - not yet configured" << std::endl;
            }
            else{
                flux[electron_index] = 0.0;
            }

            //Oxygen anions
            if(this->preceding_layer_conducts_oxygen_anions){
                
                //Extract values
                const PetscScalar mu_o_P = center_node[ion_index];
                const PetscScalar mu_o_W = west_node[ion_index];

                //Calculate conductivity
                const double sigma_o_P = this->calculate_ion_conductivity(center_node);
                //It is assumed that sigma_o_P equals sigma_o_w, which is valid as long as sigma_o_P
                //is not a function of mu_o_P or mu_e_P (the preceding layer is dense, for the concentrations
                //a Neumann boundary condition applies (slope = 0) and thus P ~ w)
                const double sigma_o_W = this->preceding_layer->calculate_ion_conductivity(west_node);
                //It is assumed that sigma_o_W equals sigma_o_w, because the conductivity is assumed to
                //change little from W to w

                //Calculate Delta_x
                const double Delta_x_Pw = this->x_center[0] - this->x_face[0];
                const double Delta_x_wW = this->preceding_layer->get_final_face_to_center_width();
                const double quotient = (sigma_o_P*Delta_x_wW)/(sigma_o_W*Delta_x_Pw);

                flux[ion_index] =   -((sigma_o_W)/(this->z[1]*this->F))*
                                                        (1.0/(1.0+quotient))*
                                                        (mu_o_W + quotient*mu_o_P - (1.0+quotient)*mu_o_W)/(Delta_x_wW);
            }
            else{
                flux[ion_index] = 0.0;
            }

        }
        else{ //Porous-porous interface
            const PetscScalar* c_i_P = center_node;
            PetscScalar c_i_w[this->number_of_gas_species], nn[this->number_of_gas_species];
            const PetscScalar* c_i_W = west_node;

            const double Delta_x_wW = this->preceding_layer->get_final_face_to_center_width();
            PetscScalar c_i_w_av[this->number_of_gas_species];
            this->calculate_concentration_face(c_i_W,c_i_w,c_i_w_av);
            PetscScalar c_i_P_av[this->number_of_gas_species];
            this->calculate_concentration_face(c_i_w,c_i_P,c_i_P_av);
 
            
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
            const double Delta_x_Pw = this->x_center[x_index] - this->x_face[x_index];
            const double quotient = (sigma_e_P*Delta_x_wW)/(sigma_e_W*Delta_x_Pw);

            flux[electron_index] =  -((sigma_e_P)/(this->z[0]*this->F))*
                                                (mu_e_P - (quotient*mu_e_P + mu_e_W)/(1.0+quotient))/(Delta_x_Pw);

            //Oxygen anions
            const PetscScalar mu_o_W = west_node[ion_index];
            const PetscScalar mu_o_P = center_node[ion_index];
            const double sigma_o_P = this->calculate_ion_conductivity(center_node);
            const double sigma_o_W = this->preceding_layer->calculate_ion_conductivity(west_node);
            const double quotient_o = (sigma_o_P*Delta_x_wW)/(sigma_o_W*Delta_x_Pw);

            flux[ion_index] =  -((sigma_o_P)/(this->z[1]*this->F))*
                                                (mu_o_P - (quotient_o*mu_o_P + mu_o_W)/(1.0+quotient_o))/(Delta_x_Pw);
        }
    }


};

void RootSubFunctionCGO::calculate_flux(    const PetscScalar* center_node, const int x_index, PetscScalar* flux,
                                            const int electron_index, const int ion_index){
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
        PetscScalar c_i_e[this->number_of_gas_species], nn[this->number_of_gas_species];
        this->calculate_concentration_face(c_i_P,c_i_bc,c_i_e);
        this->calculate_gas_flux(c_i_bc,c_i_P,c_i_e,this->x_face[x_index+1],this->x_center[x_index], nn);
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
        const double sigma_o_P = this->calculate_ion_conductivity(center_node); 
        flux[ion_index] = 0.0;// -(sigma_o_P/(this->z[1]*this->F))*
                                             //((mu_o_P- 0.0)/(this->x_face[x_index+1] - this->x_center[x_index]));


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

const double RootSubFunctionCGO::get_final_face_to_center_width() const{
    return this->x_face[this->nodes] - this->x_center[this->nodes-1];
};


void RootSubFunctionCGO::calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const{


    source[this->oxygen_index] = 0.0; //Oxygen consumption
    source[electron_index] = 0.0; //Electrons
    source[ion_index] = 0.0; //Oxygen anions
};


const double RootSubFunctionCGO::calculate_electron_conductivity(const PetscScalar* node) const{
    //?? set up storage of calculation results

    //?? implement no species situation

    // Data from Wang 2000
    const double mu_i = 1/this->temperature * std::pow(10,-3571/this->temperature +2.3667);
    const double sigma0 = mu_i*4*1.602e-19*0.2/160e-24;
    const double sigma = 1e3*sigma0*pow(1-this->porosity,1.5);
    return sigma;
};

const double RootSubFunctionCGO::calculate_ion_conductivity(const PetscScalar* node) const{

    // Data from Wang 2000
    const double mu_i = 1/this->temperature * std::pow(10,-3571/this->temperature +2.3667);
    const double sigma0 = mu_i*4*1.602e-19*0.2/160e-24;
    const double sigma = 1e2*sigma0*pow(1-this->porosity,1.5);

    return sigma;
};
