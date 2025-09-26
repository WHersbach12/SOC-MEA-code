#include "root_sub_function.h"


class RootSubFunctionCeSSZ : public RootSubFunction{
private:
    int rank;
    double temperature;
    int nodes;
    std::vector<double> x_center;
    std::vector<double> x_face;

    int preceding_layer_oxygen_anion_index;

    const double z[1] = {-2.0}; //Species charge: oxygen 

public:
    RootSubFunctionCeSSZ( const int _rank, const std::vector<double>& _x_center, const std::vector<double>& _x_face, const double _temperature,
                        const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration);

    ~RootSubFunctionCeSSZ() override{};

    void calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node,
                        const int x_index, PetscScalar* flux, const int electron_index, const int ion_index) override;

    void calculate_flux(const PetscScalar* center_node, const int x_index, PetscScalar* flux, 
                        const int electron_index, const int ion_index) override;    
    //The flux at the air compartment-electrode interface is calculated here (at x = 0)                        

    const double get_final_face_to_center_width() const override;

    void calculate_source_term(const PetscScalar* unit, PetscScalar* source, 
                                const int electron_index, const int ion_index) const override;

    const double calculate_ion_conductivity(const PetscScalar* unit) const override;

protected:


};