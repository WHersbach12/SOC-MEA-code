#include "root_sub_function_porous.h"

class RootSubFunctionMetalSupport : public RootSubFunctionPorous{
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

    const double z[1] = {-1};

public:
    RootSubFunctionMetalSupport(const int _rank, const std::vector<double>& _x_center,const std::vector<double>& _x_face,
                                const std::vector<std::string>& gas_species, const double temperature, const bool _galvanostatic,
                                const double _applied_I_or_V, const double area,
                                const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration);

    ~RootSubFunctionMetalSupport() override{};

    void calculate_flux(const PetscScalar* center_node, const PetscScalar* west_node, const int x_index, PetscScalar* flux,
                        const int electron_index, const int ion_index) override;

    void calculate_flux(const PetscScalar* center_node,  const int x_index, PetscScalar* flux, const int electron_index, const int ion_index) override; 

    const double get_final_face_to_center_width() const override;

    void calculate_source_term(const PetscScalar* node, PetscScalar* source, const int electron_index, const int ion_index) const override;

    const double calculate_electron_conductivity(const PetscScalar* unit) const override;

private:

};
