#include "root_sub_function.h"

#pragma once

class RootSubFunctionPorous : public RootSubFunction{
protected:
    const double p_atm = 101325.0; //Pa
    
    std::vector<std::string> gas_species;
    int number_of_gas_species;

    double temperature;
    double area;
    double porosity;
    double tortuosity;
    double particle_diameter;

    double* molar_mass; //kg mol^-1
    double* diffusion_volume; //-
    double* knudsen_coefficient;
    double permeability;
    double viscosity; //Pa s
    
    const GasCompartment* gas_compartment; //The gas compartment the layer has an interface with. Not set for layers 
    //that are not in contact with a gas compartment

    //Dusty Gas Model
    PetscScalar H[20][20]; //C array of H
    PetscScalar D_kl[20][20]; //C array of binary diffusion coefficients (upper triangular)
    PetscScalar D_kl_constant_quotient[20][20]; //C array of the constant quotient of the binary diffusion
    //coefficient equation - Fuller et al., 1966 (upper triangular)
  

public:
    RootSubFunctionPorous(  const std::vector<std::string>& _gas_species, const double _temperature, const double _area,
                            const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration);

    ~RootSubFunctionPorous() override;

    void set_gas_compartment(const GasCompartment* _gas_compartment) override;
    
    void set_effective_binary_diffusion_coefficients_constant_quotient(const double porosity, const double tortuosity);

    void update_effective_binary_diffusion_coefficients(const PetscScalar* c_i_w);

    //DGM methods
    void calculate_gas_flux(const PetscScalar* c_i_P, const PetscScalar* c_i_W, const PetscScalar* c_i_w,
                            const double x_center, const double x_west, PetscScalar* Flux);

    void H_Matrix_Inverse(PetscScalar H[20][20], PetscScalar D_DGM[20][20]);

protected:
    //Utility
    const double calculate_pressure(const PetscScalar* c_i) const;

    void calculate_mole_fraction(const PetscScalar* c_i, PetscScalar* x_i) const;

    void calculate_concentration_face(const PetscScalar* c_i_W, const PetscScalar* c_i_P, PetscScalar* c_i_w) const;

private:

    PetscErrorCode construct_class_petsc();

    //Internal DGM methods


};