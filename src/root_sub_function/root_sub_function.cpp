#include "root_sub_function.h"


RootSubFunction::RootSubFunction(){

};




//Dummy functions
void RootSubFunction::set_gas_compartment(const GasCompartment* _gas_compartment){
    std::cerr << "Warning: calling the dummy function set_gas_compartment()" << std::endl;
};

const double RootSubFunction::calculate_electron_conductivity(const PetscScalar* unit) const{
    std::cerr << "Warning: calling the dummy function calculate_electron_conductivity()" << std::endl;
    return 0.0;
};

const double RootSubFunction::calculate_ion_conductivity(const PetscScalar* unit) const{
    std::cerr << "Warning: calling the dummy function calculate_ion_conductivity()" << std::endl;
    return 0.0;
};

void RootSubFunction::calculate_gas_flux( const PetscScalar* c_i_P, const PetscScalar* c_i_W, const PetscScalar* c_i_w,
                                                    const double x_center, const double x_west, PetscScalar* Flux){
    std::cerr << "Warning: calling the dummy function calculate_gas_flux()" << std::endl;
};

