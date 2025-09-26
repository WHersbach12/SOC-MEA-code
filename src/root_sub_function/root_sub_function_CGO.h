#include "root_sub_function_porous.h"

class RootSubFunctionCGO : public RootSubFunctionPorous{
private:
    int rank;
    bool galvanostatic;
    double applied_I_or_V;

    int nodes;
    std::vector<double> x_center;
    std::vector<double> x_face;

    bool preceding_layer_is_dense;
    bool preceding_layer_conducts_electrons;
    bool preceding_layer_conducts_oxygen_anions;
    int preceding_layer_electron_index;
    int preceding_layer_oxygen_anion_index;

    int oxygen_index = -1;
    const double z[2] = {-1.0,-2.0}; //Species charge: MIEC electron, MIEC oxygen anion
    const double V_m = 35.5e-6; //Molar volume [m^3/mol], Matsuzaki et al. (2011)

    const double particle_area = 1.0e12;
    const double alpha = 0.5;

public:
    RootSubFunctionCGO(const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
                        const std::vector<std::string>& gas_species, const double temperature, const double area, const bool _galvanostatic,
                        const double _applied_I_or_V, const int layer_index, 
                        const std::shared_ptr<CellConfiguration>& cell_configuration);

    ~RootSubFunctionCGO() override;    

    void calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, PetscScalar* flux,
                        const int electron_index, const int ion_index) override;

    void calculate_flux(const PetscScalar* center_node,  const int x_index, PetscScalar* flux, const int electron_index, const int ion_index) override; 

    const double get_final_face_to_center_width() const override;

    void calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const override;

    const double calculate_electron_conductivity(const PetscScalar* unit) const override;

    const double calculate_ion_conductivity(const PetscScalar* unit) const override;

private:

    
};