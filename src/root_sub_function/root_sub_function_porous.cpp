#include "root_sub_function_porous.h"

RootSubFunctionPorous::RootSubFunctionPorous(   const std::vector<std::string>& _gas_species, const double _temperature, const double _area,
                                                const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration):
gas_species(_gas_species),temperature(_temperature),area(_area){

    //Derive fields from constructor input
    this->number_of_gas_species = this->gas_species.size();

    //Extract porous material properties (not stored - only need in constructor)
    // const double porosity = cell_configuration->get_layer_porosity(layer_index);
    this->porosity = cell_configuration->get_layer_porosity(layer_index);
    this->tortuosity = cell_configuration->get_layer_tortuosity(layer_index);
    this->particle_diameter = cell_configuration->get_layer_particle_diameter(layer_index);
    const double pore_diameter = (this->porosity * particle_diameter)/3/(1-this->porosity);

    //Use cell_configuration object to fill in gas species properties
    this->molar_mass = new double[this->number_of_gas_species];
    this->diffusion_volume = new double[this->number_of_gas_species];
    this->knudsen_coefficient = new double[this->number_of_gas_species];
    for(int i = 0; i < this->number_of_gas_species; i++){
        this->molar_mass[i] = cell_configuration->get_molar_mass(gas_species[i]);
        this->diffusion_volume[i] = cell_configuration->get_diffusion_volume(gas_species[i]);
        
        //Calculate the Knudsen coefficient
        const double argument = (8.0*this->R*this->temperature)/(M_PI*this->molar_mass[i]);
        this->knudsen_coefficient[i] = (1.0/3.0)*(porosity*pore_diameter/tortuosity)*std::sqrt(argument);
    }

    //Calculate the permeability (Kozeny-Carman relationship, Zhu and Kee, 2006)
    const double num = pow(porosity,3)*pow(particle_diameter,2);
    const double den = 72.0*tortuosity*pow(1.0-porosity,2);
    this->permeability = num/den;

    //Calculate the viscosity //?? implement correlations De Falco, Wilke
    this->viscosity = 43.32e-6;

    //Calculate the constant quotient of the (effective) binary diffusion coefficients equation
    // this->D_kl_constant_quotient = new double[this->number_of_gas_species][this->number_of_gas_species];
    this->set_effective_binary_diffusion_coefficients_constant_quotient(porosity,tortuosity);

    //?? check single species case!

    //Perform class construction specific to PetSc
    const int petsc_error = construct_class_petsc();

};


RootSubFunctionPorous::~RootSubFunctionPorous(){
    
    //Clean-up C arrays
    delete[] this->molar_mass;
    delete[] this->diffusion_volume;
    delete[] this->knudsen_coefficient;

    // delete[] this->D_kl;    
    // delete[] this->D_kl_constant_quotient;

};


PetscErrorCode RootSubFunctionPorous::construct_class_petsc(){

    //Create matrix (DGM)
    // this->H = new PetscScalar[this->number_of_gas_species][this->number_of_gas_species]; //Logically 2D, row-oriented
    // this->D_kl = new PetscScalar[this->number_of_gas_species][this->number_of_gas_species]; //Logically 2D, row-oriented

    PetscFunctionReturn(0);
};


void RootSubFunctionPorous::set_gas_compartment(const GasCompartment* _gas_compartment){
    this->gas_compartment = _gas_compartment;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// DGM
void RootSubFunctionPorous::calculate_gas_flux( const PetscScalar* c_i_P, const PetscScalar* c_i_W, const PetscScalar* c_i_w,
                                                const double x_center, const double x_west, PetscScalar* Flux){
    // Fill the H matrix
    PetscScalar x_w[this->number_of_gas_species], x_W[this->number_of_gas_species], x_P[this->number_of_gas_species];
    this->calculate_mole_fraction(c_i_w, x_w);
    this->update_effective_binary_diffusion_coefficients(c_i_w);

    for(int k = 0; k < this->number_of_gas_species; k++){ //Iterate over all rows 
        for(int l = 0; l < this->number_of_gas_species; l++){ //Iterate over columns (upper triangular matrix)
            if(k == l){
                double x_over_Dkl = 0;
                for(int j = 0; j < this->number_of_gas_species; j++){
                    if(j != k){
                        x_over_Dkl += x_w[j]/this->D_kl[k][j];
                    }
                }
                H[k][l] = 1/this->knudsen_coefficient[k] + x_over_Dkl;//sum of x_k/D_kl
            }
            else{
                H[k][l] = x_w[l]/this->D_kl[k][l];
            }
        }
    }

    PetscScalar D_DGM[20][20];
    this->H_Matrix_Inverse(H,D_DGM);

    // Calculate the fluxes
    const double P_P = this->calculate_pressure(c_i_P);
    const double P_W = this->calculate_pressure(c_i_W);

    for(int k = 0; k < this->number_of_gas_species; k++){
        PetscScalar Diff_flux = 0;
        PetscScalar DGM_conv = 0;
        for(int l = 0; l < this->number_of_gas_species; l++){
            Diff_flux += D_DGM[k][l]*(c_i_P[l] - c_i_W[l])/(x_center - x_west);
            DGM_conv += D_DGM[k][l]*c_i_w[l]/this->knudsen_coefficient[l];
            
        }
        Flux[k] = - Diff_flux - DGM_conv*this->permeability/this->viscosity * (P_P - P_W)/(x_center - x_west);
    }
};

void RootSubFunctionPorous::H_Matrix_Inverse(PetscScalar H[20][20], PetscScalar D_DGM[20][20]){
    // Matrix is inverted by doing LU factorization. This creates two triangular matrices: an upper and a lower.
    PetscScalar L[this->number_of_gas_species][this->number_of_gas_species];
    PetscScalar U[this->number_of_gas_species][this->number_of_gas_species];
    PetscScalar L_inv[this->number_of_gas_species][this->number_of_gas_species];
    PetscScalar U_inv[this->number_of_gas_species][this->number_of_gas_species];
    PetscScalar D_DGM2[this->number_of_gas_species][this->number_of_gas_species];

    //Make the triangular matrices.
    for(int j = 0; j < this->number_of_gas_species; j++){
        for(int i = 0; i < this->number_of_gas_species; i++){
            L[i][j] = 0.0; U[i][j] = 0.0; L_inv[i][j] = 0.0; U_inv[i][j] = 0.0;
            if(j==0){L[i][j] = H[i][j];}
            if(i==0){U[i][j] = H[i][j]/H[0][0];}
            else if(i==j){U[i][j] = 1;}
        }
    }

    for(int j = 1; j < this->number_of_gas_species; j++){
        for(int i = j; i < this->number_of_gas_species; i++){
            L[i][j] = H[i][j];
            for(int k = 0; k < j; k++){L[i][j] -= L[i][k]*U[k][j];}
        }
        for(int k = j+1; k < this->number_of_gas_species; k++){
            U[j][k] = H[j][k]/L[j][j];
            for(int i = 0; i < j; i++){U[j][k] -= L[j][i]*U[i][k]/L[j][j];}
        }
    }

    //Compute the inverses of L
    for(int i = 0; i < this->number_of_gas_species; i++){
        L_inv[i][i] = 1/L[i][i];
        for(int j = 0; j < i; j++){
            for(int k = j; k < i; k++){
                L_inv[i][j] += L[i][k]*L_inv[k][j];
            }
            L_inv[i][j] = -L_inv[i][j]/L[i][i];
        }
    }
    //Compute the inverses of U
    for(int i = 0; i < this->number_of_gas_species; i++){
        U_inv[i][i] = 1/U[i][i];
        for(int j = 0; j < i; j++){
            for(int k = j; k < i; k++){
                U_inv[j][i] += U[k][i]*U_inv[j][k];
            }
            U_inv[j][i] = -U_inv[j][i]/U[i][i];
        }
    }
    // Calculate the inverse H matrix by LU factorization.
    for(int i = 0; i < this->number_of_gas_species; i++){
        for(int j = 0; j < this->number_of_gas_species; j++){
            D_DGM[i][j] = 0.0;
            for(int k = 0; k < this->number_of_gas_species; k++){
                D_DGM[i][j] += U_inv[i][k]*L_inv[k][j];
            }
        }
    }
};


void RootSubFunctionPorous::set_effective_binary_diffusion_coefficients_constant_quotient(const double porosity, const double tortuosity){
    //constant quotient: (epsilon/tau)*sqrt(1/M_i + 1/M_j)/((V_i^1/3 + V_j^1/3)^2)

    const double prefactor = porosity/tortuosity;
    for(int k = 0; k < this->number_of_gas_species; k++){ //Iterate over all rows 
        for(int l = 0; l < this->number_of_gas_species; l++){ //Iterate over columns (upper triangular matrix)
           
            const double num = std::sqrt(1e-3/this->molar_mass[k] + 1e-3/this->molar_mass[l]);
            const double den = std::pow(this->diffusion_volume[k],1.0/3.0) + std::pow(this->diffusion_volume[l],1.0/3.0);
            const double den2 = den*den;

            this->D_kl_constant_quotient[k][l] = prefactor*num/den2;
        }
    }

};

void RootSubFunctionPorous::update_effective_binary_diffusion_coefficients(const PetscScalar* c_i_w){

    const double p = this->calculate_pressure(c_i_w); //Pressure in Pa
    const double prefactor = (1.0e-7*std::pow(this->temperature,1.75))/(p/this->p_atm);

    for(int k = 0; k < this->number_of_gas_species; k++){ //Iterate over all rows 
        for(int l = 0; l < this->number_of_gas_species; l++){ //Iterate over columns (upper triangular matrix)
            this->D_kl[k][l] = prefactor*this->D_kl_constant_quotient[k][l];
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Utility

const double RootSubFunctionPorous::calculate_pressure(const PetscScalar* c_i) const{
    double c = 0.0; //Total concentration
    for(int i = 0; i < this->number_of_gas_species; i++){c += c_i[i];}
    return c*this->R*this->temperature;
};

void RootSubFunctionPorous::calculate_mole_fraction(const PetscScalar* c_i, PetscScalar* x_i) const{
    double c = 0.0; //Total concentration
    for(int i = 0; i < this->number_of_gas_species; i++){c += c_i[i];}
    for(int i = 0; i < this->number_of_gas_species; i++){x_i[i] = c_i[i]/c;}
};

void RootSubFunctionPorous::calculate_concentration_face(const PetscScalar* c_i_W, const PetscScalar* c_i_P, PetscScalar* c_i_w) const{
    for(int i = 0; i < this->number_of_gas_species; i++){c_i_w[i] = (c_i_W[i] + c_i_P[i])/2.0;}
};

