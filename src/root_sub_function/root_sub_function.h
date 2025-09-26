#include "../cell_configuration.h"

#pragma once

class RootSubFunction{
protected:
    const double F = 96485.33; //Faraday's constant [C/mol]
    const double R = 8.3145; //Gas constant [J/mol]

    std::shared_ptr<RootSubFunction> preceding_layer;

public:
    RootSubFunction();

    virtual ~RootSubFunction(){};

    inline void set_preceding_layer(const std::shared_ptr<RootSubFunction>& _preceding_layer){
        this->preceding_layer = _preceding_layer;
    }

    virtual void set_gas_compartment(const GasCompartment* _gas_compartment);

    //Note that all flux are calculated in the positve x-direction
    virtual void calculate_flux(const PetscScalar* unit_center, const PetscScalar* unit_west, const int x_index, 
                                PetscScalar* flux, const int electron_index, const int ion_index) = 0;

    virtual void calculate_flux(const PetscScalar* unit_center, const int x_index, PetscScalar* flux, 
                                const int electron_index, const int ion_index) = 0;

    virtual const double get_final_face_to_center_width() const = 0;
    //Used for boundary conditions

    virtual void calculate_source_term( const PetscScalar* unit, PetscScalar* source, 
                                        const int electron_index, const int ion_index) const = 0;
    //All source terms

    //Virtual functions that can, but do not have to, be overriden
    virtual const double calculate_electron_conductivity(const PetscScalar* unit) const;

    virtual const double calculate_ion_conductivity(const PetscScalar* unit) const;
    //DGM method

    virtual void calculate_gas_flux( const PetscScalar* c_i_P, const PetscScalar* c_i_W, const PetscScalar* c_i_w,
                                            const double x_center, const double x_west, PetscScalar* Flux); 

private:


};