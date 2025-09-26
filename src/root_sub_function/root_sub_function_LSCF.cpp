#include "root_sub_function_LSCF.h"

RootSubFunctionLSCF::RootSubFunctionLSCF(   const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
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
        std::cerr << "Warning: oxygen species not found for LSCF - not compatible with electron conductivity correlation (RootSubFunctionLSCF)" << std::endl;
    }

    //Percolation probability based on [Bertei 2011] A comparative study and an extended theory of percolation for random packings of rigid spheres
    // and [Bertei 2011] Percolation theory in SOFC composite electrodes
    const double nr_particles = (1-this->porosity)/(3.1416/3*this->particle_diameter*this->particle_diameter);
    this->particle_area = 6*(1-this->porosity)/this->particle_diameter;//3.1416*nr_particles*this->particle_diameter*this->particle_diameter;

    const double dH_O2 = 0 + 31.32234*this->temperature/1000 - 20.23531*pow(this->temperature/1000,2)/2 + 57.8664*pow(this->temperature/1000,3)/3 - 
                        36.50624*pow(this->temperature/1000,4)/4 - 0.007374*pow(this->temperature/1000,-1);
    const double dS_O2 = 205.152;
    this->mu_O2_0 = dH_O2 - this->temperature*dS_O2;

    this->set_effective_binary_diffusion_coefficients_constant_quotient(porosity,tortuosity);
};

RootSubFunctionLSCF::~RootSubFunctionLSCF(){};    

void RootSubFunctionLSCF::calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, PetscScalar* flux, 
                    const int electron_index, const int ion_index){
    //Calculate the flux at the western cell face, an array with length equal to center_node_length
    
    //x_index: x index of cell center P

    if(x_index > 0){ //Flux at all cell faces inside LSCF
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
    else{ //As this version of the overload is called, there must exist a preceding node and LSCF does not interface the air compartment
        
        //?? Not configured yet
        std::cerr << "Attempt to calculate the flux at the western interface between LSCF and a preceding layer - not yet configured" << std::endl;

    }


};

void RootSubFunctionLSCF::calculate_flux(const PetscScalar* center_node, const int x_index, PetscScalar* flux, const int electron_index, const int ion_index){
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
        this->calculate_gas_flux(c_i_P,c_i_bc,c_i_w,this->x_center[0],this->x_face[0], nn);
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
                                                ((mu_e_P - E_cell)/(this->x_center[0] - this->x_face[0]));
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

const double RootSubFunctionLSCF::get_final_face_to_center_width() const{
    return this->x_face[this->nodes] - this->x_center[this->nodes-1];
};


void RootSubFunctionLSCF::calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const{

    //?? implement no species situation

    //Calculate the surface overpotential
    const double mu_e = node[electron_index];
    const double mu_o = node[ion_index];
    const double p_O2 = node[this->oxygen_index]*this->R*this->temperature/this->p_atm; 
    const double Delta_chi = (-2*mu_e + mu_o + 0.5*this->mu_O2_0 + this->R*this->temperature*std::log(std::sqrt(p_O2)))/(2*this->F);
    
    const double i0 = 1994.64*std::exp(-149861/this->R/this->temperature)*this->particle_area*2*this->F*std::pow(p_O2,0.25);

    //Calculate the Butler-Volmer equation
    const double i = i0 * (std::exp(this->alpha*2*this->F*Delta_chi/(this->R*this->temperature)) - 
                    std::exp(-(1.0-this->alpha)*2*this->F*Delta_chi/(this->R*this->temperature)) );

    //?? temporary
    source[this->oxygen_index] = -i/(4*this->F); //Oxygen consumption
    source[electron_index] = -i; //Electrons
    source[ion_index] = i; //Oxygen anions
};


const double RootSubFunctionLSCF::calculate_electron_conductivity(const PetscScalar* node) const{
    //?? set up storage of calculation results

    //?? implement no species situation

    const double p_O2 = node[this->oxygen_index]*this->R*this->temperature/this->p_atm; 
    // const double p_O2 = 0.21;
    const double log_p_O2 = std::log(p_O2)/std::log(10.0);
    const double log_sigma = -0.0237*log_p_O2*log_p_O2 + 0.00034*log_p_O2 + 4.8126; //Bouwmeester et al, 1073.15 K
    const double sigma0 = std::pow(10.0,log_sigma);
    // const double sigma = sigma0*1999*pow(this->particle_diameter*1e6,3.305)*pow(1-this->porosity,1.5);
    const double sigma = sigma0*pow(1-this->porosity,1.5);
    return sigma;
};

const double RootSubFunctionLSCF::calculate_ion_conductivity(const PetscScalar* node) const{
    const double p_O2 = node[this->oxygen_index]*this->R*this->temperature/this->p_atm; 
    // const double p_O2 = 0.21;
    const double log_p_O2 = std::log(p_O2)/std::log(10.0);
    const double log_D = -0.1765*log_p_O2*log_p_O2 -0.2724*log_p_O2 - 9.2256; //Bouwmeester et al, 1073.15 K
    const double D = std::pow(10.0,log_D);
    const double ddelta_dlnpO2 = -3.36260e-5*this->temperature + 2.59403e-2; //Bouwmeester et al.
    const double sigma0 = -D*((8.0*this->F*this->F)/(this->R*this->temperature*this->V_m))*ddelta_dlnpO2;
    // const double sigma = sigma0*1999*pow(this->particle_diameter*1e6,3.305)*pow(1-this->porosity,1.5);
    const double sigma = sigma0*pow(1-this->porosity,1.5);

    return sigma;
};
