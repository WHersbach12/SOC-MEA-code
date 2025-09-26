#include "root_sub_function_LSMYSZ.h"

RootSubFunctionLSMYSZ::RootSubFunctionLSMYSZ(   const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
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

    //Percolation probability based on [Bertei 2011] A comparative study and an extended theory of percolation for random packings of rigid spheres
    // and [Bertei 2011] Percolation theory in SOFC composite electrodes
    const double ratio_dp = 1.0;
    const double S_LSM = this->Phi_LSM/(this->Phi_LSM + ratio_dp*(1 - this->Phi_LSM));
    const double S_YSZ = this->Phi_YSZ/(this->Phi_YSZ + ratio_dp*(1 - this->Phi_YSZ));
    this->Percolation_LSM = 1 - pow((4.236 - 6*S_LSM)/(2.472) , 3.7);
    this->Percolation_YSZ = 1 - pow((4.236 - 6*S_YSZ)/(2.472) , 3.7);
    const double zeta_LSM = this->Phi_LSM/(this->Phi_LSM + pow(ratio_dp,3.0)*(1 - this->Phi_LSM));
    const double zeta_YSZ = this->Phi_YSZ/(this->Phi_YSZ + pow(ratio_dp,3.0)*(1 - this->Phi_YSZ));

    const double nr_particles = (1-this->porosity)/(3.1416/6*(zeta_LSM*pow(this->particle_diameter,3)+zeta_YSZ*pow(this->particle_diameter,3)));
    const double r_c = this->particle_diameter/2*sin(15*3.1416/180);
    const double ZLSM_YSZ = 2*(3+(6-3)*pow(this->particle_diameter,2)/((zeta_LSM+zeta_YSZ)*pow(this->particle_diameter,2)));
    //const double ZNi_YSZ = (3*(2-sqrt(3))*6*(ratio_dp+1))/(1 + ratio_dp - sqrt(ratio_dp*(ratio_dp+2)));

    this->TPB_length = nr_particles*zeta_LSM*zeta_YSZ*Percolation_LSM*Percolation_YSZ*ZLSM_YSZ*2*3.1416*r_c;

    const double dH_O2 = 0 + 31.32234*this->temperature/1000 - 20.23531*pow(this->temperature/1000,2)/2 + 57.8664*pow(this->temperature/1000,3)/3 - 
                        36.50624*pow(this->temperature/1000,4)/4 - 0.007374*pow(this->temperature/1000,-1);
    const double dS_O2 = 205.152;
    this->mu_O2_0 = dH_O2 - this->temperature*dS_O2;


    this->set_effective_binary_diffusion_coefficients_constant_quotient(this->porosity,this->tortuosity);
};

RootSubFunctionLSMYSZ::~RootSubFunctionLSMYSZ(){};    

void RootSubFunctionLSMYSZ::calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, PetscScalar* flux, 
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
            const PetscScalar mu_o_W = west_node[ion_index];
            const PetscScalar mu_o_P = center_node[ion_index];
            const double sigma_o_P = this->calculate_ion_conductivity(center_node);
            const double sigma_o_W = this->preceding_layer->calculate_ion_conductivity(west_node);
            const double quotient_o = (sigma_o_P*Delta_x_wW)/(sigma_o_W*Delta_x_Pw);

            flux[ion_index] =  -((sigma_o_P)/(this->z[1]*this->F))*
                                                (mu_o_P - (quotient_o*mu_o_P + mu_o_W)/(1.0+quotient_o))/(Delta_x_Pw);           
        }
        else{ //Porous - porous interfaces
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
            if(preceding_layer_conducts_oxygen_anions){
                //Oxygen anions
                const PetscScalar mu_o_W = west_node[ion_index];
                const PetscScalar mu_o_P = center_node[ion_index];
                const double sigma_o_P = this->calculate_ion_conductivity(center_node);
                const double sigma_o_W = this->preceding_layer->calculate_ion_conductivity(west_node);
                const double quotient_o = (sigma_o_P*Delta_x_wW)/(sigma_o_W*Delta_x_Pw);

                flux[ion_index] =  -((sigma_o_P)/(this->z[1]*this->F))*
                                                    (mu_o_P - (quotient_o*mu_o_P + mu_o_W)/(1.0+quotient_o))/(Delta_x_Pw);
            }
            else{
                flux[ion_index] = 0.0;
            }
        }
    }


};

void RootSubFunctionLSMYSZ::calculate_flux(const PetscScalar* center_node, const int x_index, PetscScalar* flux, const int electron_index, const int ion_index){
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

const double RootSubFunctionLSMYSZ::get_final_face_to_center_width() const{
    return this->x_face[this->nodes] - this->x_center[this->nodes-1];
};

void RootSubFunctionLSMYSZ::calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const{

    //?? implement no species situation
    double x_i[this->number_of_gas_species];
    this->calculate_mole_fraction(node,x_i);
    //Calculate the surface overpotential
    const double mu_e = node[electron_index];
    const double mu_o = node[ion_index];
    const double p_O2 = node[this->oxygen_index]*this->R*this->temperature/this->p_atm; 
    const double p_O2_eq = this->gas_compartment->inlet_mole_fractions[this->oxygen_index];

    const double Delta_chi = (-2*mu_e + mu_o + 0.5*this->mu_O2_0 + 0.5*this->R*this->temperature*std::log(x_i[this->oxygen_index]))/(2*this->F);
    
    //Kinetic data 
    const double i0 = 1e-9*std::exp(-10e3/this->R/this->temperature)*this->TPB_length*2*this->F*std::pow(p_O2,0.25);
    //Calculate the Butler-Volmer equation
    // const double i = 10.0*this->F;
    const double i = i0 * (std::exp(this->alpha*2*this->F*Delta_chi/(this->R*this->temperature)) - 
                    std::exp(-(1.0-this->alpha)*2*this->F*Delta_chi/(this->R*this->temperature)) );

    //?? temporary
    source[this->oxygen_index] = -i/(4*this->F); //Oxygen consumption
    source[electron_index] = -i; //Electrons
    source[ion_index] = i; //Oxygen anions
};


const double RootSubFunctionLSMYSZ::calculate_electron_conductivity(const PetscScalar* node) const{

    // [Zhu 2008] Modelling distributed charge-transfer proceses in SOFC membrane electrode assemblies
    double sigma0 = 8.855e5/this->temperature * exp(-1082.5/this->temperature);
    double sigma = 1e2*sigma0*pow((1 - this->porosity)*this->Phi_LSM*this->Percolation_LSM, 1.5);

    return sigma;
};

const double RootSubFunctionLSMYSZ::calculate_ion_conductivity(const PetscScalar* node) const{

    // [Zhu 2008] Modelling distributed charge-transfer proceses in SOFC membrane electrode assemblies
    // double sigma0 = 3.34e2/this->temperature * exp(-10300/this->temperature);
    double sigma0 = 3.6e5/this->temperature * exp(-80.116e3/this->R/this->temperature);
    double sigma = 1e2*sigma0*pow((1 - this->porosity)*this->Phi_YSZ*this->Percolation_YSZ, 1.5);

    return sigma;
};
