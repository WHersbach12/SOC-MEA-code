#include "root_sub_function_porous.h"

class RootSubFunctionNiYSZ : public RootSubFunctionPorous{
private:
    int rank;

    bool galvanostatic;
    double applied_I_or_V;
    int nodes;
    std::vector<double> x_center;
    std::vector<double> x_face;

    int hydrogen_index = -1;
    int water_index = -1;
    
    bool preceding_layer_is_dense;
    bool preceding_layer_conducts_electrons;
    bool preceding_layer_conducts_oxygen_anions;
    int preceding_layer_electron_index;
    int preceding_layer_oxygen_anion_index;

    const double z[2] = {-1.0,-2.0}; //Species charge: pure-conductor electron, pure-conductor oxygen anion

    double Percolation_Ni;
    double Percolation_YSZ;
    const double Phi_YSZ = 0.54;
    const double Phi_Ni = 1-Phi_YSZ;
    const double alpha = 0.5;
    double TPB_length;
    double mu_H2_0;
    double mu_H2O_0;

public:
    RootSubFunctionNiYSZ(   const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
                            const std::vector<std::string>& gas_species, const double temperature, const bool _galvanostatic,
                            const double _applied_I_or_V, const double area,
                            const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration);

    ~RootSubFunctionNiYSZ() override{};

    void calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, PetscScalar* flux,
                        const int electron_index, const int ion_index) override;

    void calculate_flux(const PetscScalar* center_node,  const int x_index, PetscScalar* flux, const int electron_index, const int ion_index) override; 

    const double get_final_face_to_center_width() const override;

    void calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const override;

    const double calculate_electron_conductivity(const PetscScalar* node) const override;

    const double calculate_ion_conductivity(const PetscScalar* node) const override;

private:

};