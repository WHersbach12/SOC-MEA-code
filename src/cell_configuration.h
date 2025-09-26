#include "dependencies.h"

#pragma once
typedef std::unordered_map<std::string,std::pair<int,int>> umap_sii; //key + (value + value)
typedef std::unordered_map<std::string,double> umap_sd; //key + value

struct Layer{
    std::string name;
    bool is_dense;  //Is the layer dense (1) or porous (0)?
    bool is_air_electrode;  //Is the layer in the air (1) or fuel (0) electrode?
    double thickness;
    int nodes;
    double porosity;
    double tortuosity;
    double particle_diameter;
    std::vector<std::string> gas_species;
    int number_of_gas_species;
    std::pair<int,int> charge_carriers; //electron, oxygen anions - main conductor first
    int number_of_charge_carriers;

    std::vector<double> x_center;
    std::vector<double> x_face;

    std::vector<int> mapping_to_next_layer; //The mapping from the dependent variables of
    //the current layer to the dependent variables of the next layer


    //Configure serialization
    template<class Archive> void serialize(Archive& ar, const unsigned int version){
        ar & name;
        ar & is_dense;
        ar & is_air_electrode;
        ar & thickness;
        ar & nodes;
        ar & porosity;
        ar & tortuosity;
        ar & particle_diameter;
        ar & gas_species;
        ar & charge_carriers;
        ar & number_of_gas_species;
        ar & number_of_charge_carriers;
        ar & x_center;
        ar & x_face;
        ar & mapping_to_next_layer;

    }

};

struct GasCompartment{
    double inlet_pressure; //Pa
    double inlet_flow_rate; //m_g^3/s
    int number_of_gas_species;
    std::vector<std::string> gas_species;
    std::vector<double> inlet_mole_fractions;

    //Configure serialization
    template<class Archive> void serialize(Archive& ar, const unsigned int version){
        ar & inlet_pressure;
        ar & inlet_flow_rate;
        ar & number_of_gas_species;
        ar & gas_species;
        ar & inlet_mole_fractions;
    }

};


//Stateful
class CellConfiguration
{
private:
    
    int layer_count;
    std::vector<Layer> layers;
    std::pair<GasCompartment,GasCompartment> gas_compartments;

    std::vector<std::string> valid_layer_names;
    umap_sii layer_charge_carriers; //Layer name + (number of electron conservation equations + number of ion conservation equations)

    std::vector<int> air_electrode_layer_indices; //Indices of layers belonging to the air electrode (before the dense section)
    std::vector<int> fuel_electrode_layer_indices; //Indices of layers belonging to the fuel electrode (after the dense section)

    //Gas properties
    std::vector<std::string> valid_gas_species;
    umap_sd molar_mass;
    umap_sd diffusion_volume;
    
    //Configure serialization
    friend class boost::serialization::access;

    template<class Archive> void serialize(Archive& ar, const unsigned int version){
        ar & layer_count;
        ar & valid_layer_names;
        ar & layer_charge_carriers;
        ar & air_electrode_layer_indices;
        ar & fuel_electrode_layer_indices;
        ar & layers;
        ar & gas_compartments;
        ar & molar_mass;
        ar & diffusion_volume;
        ar & ion_index;
        ar & electron_index;
    }
    
public:
    int ion_index;
    int electron_index;
    std::vector<int> layer_start_index;
    
    CellConfiguration();

    void construct_class();

    void add_layer( const std::string& name, const double thickness, const int nodes); //Overload for dense layers

    void add_layer( const std::string& name, const double thickness, const int nodes, const double porosity, 
                    const double tortuosity, const double particle_diameter, const std::vector<std::string>& gas_species); //Overload for porous layers

    void add_gas_compartment(   const bool is_air, const double inlet_pressure, const double inlet_flow_rate, 
                                const std::vector<std::string>& gas_species, const std::vector<double>& inlet_mole_fractions);

    void print_configuration() const;

    void check_compatibility_gas_interfaces();

    //Calculate the length for every unit of the solution vector, the starting index for every layer and the total solution vector length
    void calculate_solution_vector_unit_lengths(std::vector<int>& node_layer_indices, int& total, 
                                                int& number_of_variables);

    //Layer
    inline const int get_layer_count() const{
        return this->layer_count;
    }

    inline const std::string& get_layer_name(const int layer) const{
        return this->layers[layer].name;
    }

    inline const bool get_layer_is_dense(const int layer) const{
        return this->layers[layer].is_dense;
    }

    inline const std::vector<std::string>& get_layer_gas_species(const int layer) const{
        return this->layers[layer].gas_species;
    }

    inline const int get_layer_nodes(const int layer) const{
        return this->layers[layer].nodes;
    }

    inline const double get_layer_porosity(const int layer) const{
        double output = -1.0;
        if(this->layers[layer].is_dense){
            std::cerr << "Warning: attempting to retrieve the porosity of a dense layer! (layer index: " << layer << ") (get_layer_porosity)" << std::endl;
        }
        else{
            output = this->layers[layer].porosity;
        }
        return output;
    }

    inline const double get_layer_tortuosity(const int layer) const{
        double output = -1.0;
        if(this->layers[layer].is_dense){
            std::cerr << "Warning: attempting to retrieve the tortuosity of a dense layer! (layer index: " << layer << ") (get_layer_tortuosity)" << std::endl;
        }
        else{
            output = this->layers[layer].tortuosity;
        }
        return output;
    }

    inline const double get_layer_particle_diameter(const int layer) const{
        double output = -1.0;
        if(this->layers[layer].is_dense){
            std::cerr << "Warning: attempting to retrieve the pore diameter of a dense layer! (layer index: " << layer << ") (get_particle_diameter)" << std::endl;
        }
        else{
            output = this->layers[layer].particle_diameter;
        }
        return output;
    }

    inline const std::pair<int,int> get_layer_charge_carriers(const int layer) const{
        return this->layers[layer].charge_carriers;
    }    

    inline const int get_layer_number_of_gas_species(const int layer) const{
        return this->layers[layer].number_of_gas_species;
    }

    inline const int get_layer_number_of_charge_carriers(const int layer) const{
        return this->layers[layer].number_of_charge_carriers;
    }

    inline const std::vector<double>& get_layer_x_center(const int layer) const{
        return this->layers[layer].x_center;
    }

    inline const std::vector<double>& get_layer_x_face(const int layer) const{
        return this->layers[layer].x_face;
    }

    inline const bool get_layer_is_air_electrode(const int layer) const{
        return this->layers[layer].is_air_electrode;
    }

    //Gas compartment
    inline const double get_gas_compartment_inlet_pressure(const bool is_air) const{
        return is_air ? this->gas_compartments.first.inlet_pressure : this->gas_compartments.second.inlet_pressure;
    }

    inline const double get_gas_compartment_inlet_flow_rate(const bool is_air) const{
        return is_air ? this->gas_compartments.first.inlet_flow_rate : this->gas_compartments.second.inlet_flow_rate;
    }

    inline const int get_gas_compartment_number_of_gas_species(const bool is_air) const{
        return is_air ? this->gas_compartments.first.number_of_gas_species : this->gas_compartments.second.number_of_gas_species;
    }

    inline const std::vector<std::string>& get_gas_compartment_gas_species(const bool is_air) const{
        return is_air ? this->gas_compartments.first.gas_species : this->gas_compartments.second.gas_species;
    }

    inline const std::vector<double>& get_gas_compartment_inlet_mole_fractions(const bool is_air) const{
        return is_air ? this->gas_compartments.first.inlet_mole_fractions : this->gas_compartments.second.inlet_mole_fractions;
    }

    inline const GasCompartment* get_gas_compartment(const bool is_air) const{
        return is_air ? &this->gas_compartments.first : &this->gas_compartments.second;
    }

    //Gas species properties
    inline const double get_molar_mass(const std::string& species_name) const{
        return this->molar_mass.at(species_name);
    }

    inline const double get_diffusion_volume(const std::string& species_name) const{
        return this->diffusion_volume.at(species_name);
    }


private:

    void generate_grid(const double x_begin, const int nodes, const double thickness, std::vector<double>& center, std::vector<double>& face) const;

};