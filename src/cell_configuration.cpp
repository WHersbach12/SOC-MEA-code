#include "cell_configuration.h"

CellConfiguration::CellConfiguration(){};

void CellConfiguration::construct_class(){
    this->layer_count = 0; //Set layer count to 0

    //Set valid layer names
    this->valid_layer_names = {
        "LSCF",
        "LSCF-CGO",
        "LSM-YSZ",
        "YSZ",
        "3YSZ",
        "Ni-YSZ",
        "Ni-CGO",
        "Support",
        "CGO",
        "LSM-YSZ",
        "LSM",
        "CeSSZ"
    };

    //Set layer_charge_carriers
    this->layer_charge_carriers.insert(std::make_pair("LSCF",std::make_pair(1,1)));
    this->layer_charge_carriers.insert(std::make_pair("LSCF-CGO",std::make_pair(1,2)));
    this->layer_charge_carriers.insert(std::make_pair("LSM-YSZ",std::make_pair(1,1)));
    this->layer_charge_carriers.insert(std::make_pair("YSZ",std::make_pair(0,1)));
    this->layer_charge_carriers.insert(std::make_pair("3YSZ",std::make_pair(0,1)));
    this->layer_charge_carriers.insert(std::make_pair("CeSSZ",std::make_pair(0,1)));
    this->layer_charge_carriers.insert(std::make_pair("Ni-YSZ",std::make_pair(1,1)));
    this->layer_charge_carriers.insert(std::make_pair("Ni-CGO",std::make_pair(1,1)));
    this->layer_charge_carriers.insert(std::make_pair("Support",std::make_pair(1,0)));
    this->layer_charge_carriers.insert(std::make_pair("CGO",std::make_pair(1,1)));
    this->layer_charge_carriers.insert(std::make_pair("LSM",std::make_pair(1,0)));

    //Check that layer_charge_carriers is set for all valid layer names
    for(int i = 0; i < this->valid_layer_names.size(); i++){
        if(this->layer_charge_carriers.count(this->valid_layer_names[i]) < 1){
            std::cerr << "Error: for layer " << this->valid_layer_names[i] << " no entry was found in layer_charge_carriers (CellConfiguration)" << std::endl;
        }
    }    


    this->valid_gas_species = {
        "N2",
        "O2",
        "H2O",
        "H2",
        "CO",
        "CO2"
    };

    //Set the molar mass in kg mol^-1    
    this->molar_mass.insert(std::make_pair("N2",28.0134e-3));
    this->molar_mass.insert(std::make_pair("O2",32.0e-3));
    this->molar_mass.insert(std::make_pair("H2O",18.01528e-3));
    this->molar_mass.insert(std::make_pair("H2",2.016e-3));
    this->molar_mass.insert(std::make_pair("CO",28.01e-3));
    this->molar_mass.insert(std::make_pair("CO2",44.01e-3));

    //Set the diffusion volumes - Fuller et al. 1966
    this->diffusion_volume.insert(std::make_pair("N2",17.9));
    this->diffusion_volume.insert(std::make_pair("O2",16.6));
    this->diffusion_volume.insert(std::make_pair("H2O",12.7));
    this->diffusion_volume.insert(std::make_pair("H2",7.07));
    this->diffusion_volume.insert(std::make_pair("CO",18.9));
    this->diffusion_volume.insert(std::make_pair("CO2",26.9));

    //Check that molar_mass and diffusion_volume is set for all valid layer names
    for(int i = 0; i < this->valid_gas_species.size(); i++){
        if(this->molar_mass.count(this->valid_gas_species[i]) < 1){
            std::cerr << "Error: for species " << this->valid_gas_species[i] << " no entry was found in molar_mass (CellConfiguration)" << std::endl;
        }
        if(this->diffusion_volume.count(this->valid_gas_species[i]) < 1){
            std::cerr << "Error: for species " << this->valid_gas_species[i] << " no entry was found in diffusion_volume (CellConfiguration)" << std::endl;
        }
    }    

}

//Overload for dense layers
void CellConfiguration::add_layer( const std::string& name, const double thickness, const int nodes){ 
    Layer new_layer;

    //Check if the layer name is valid
    bool is_valid_layer_name = 0;
    for(int i = 0; i < this->valid_layer_names.size(); i++){
        if(this->valid_layer_names[i] == name){
            is_valid_layer_name = 1;
            new_layer.name = name;
            break;
        }
    }
    if(!is_valid_layer_name){
        std::cerr << "Error: provided layer name (" << name << ") is invalid (add_layer)" << std::endl;
    }

    new_layer.is_dense = 1;
    new_layer.number_of_gas_species = 0;
    new_layer.thickness = thickness;
    new_layer.nodes = nodes;

    new_layer.charge_carriers = this->layer_charge_carriers.at(name);
    new_layer.number_of_charge_carriers = new_layer.charge_carriers.first + new_layer.charge_carriers.second;
    
    //Generate grid
    std::vector<double> center(nodes);
    std::vector<double> face(nodes+1);
    if(this->layer_count == 0){
        this->generate_grid(0.0,nodes,thickness,center,face);
        new_layer.x_center = center;
        new_layer.x_face = face;        
    }
    else{
        const int nodes_previous = this->layers[this->layer_count-1].nodes;
        const double x_begin = this->layers[this->layer_count-1].x_face[nodes_previous]; //Last face values of previous layer
        this->generate_grid(x_begin,nodes,thickness,center,face);
        new_layer.x_center = center;
        new_layer.x_face = face;        
    }

    this->layers.push_back(new_layer);
    this->layer_count++;    
}; 


//Overload for porous layers
void CellConfiguration::add_layer(  const std::string& name, const double thickness, const int nodes, const double porosity, 
                                    const double tortuosity, const double particle_diameter, const std::vector<std::string>& gas_species){ 
    Layer new_layer;

    //Check if the layer name is valid
    bool is_valid_layer_name = 0;
    for(int i = 0; i < this->valid_layer_names.size(); i++){
        if(this->valid_layer_names[i] == name){
            is_valid_layer_name = 1;
            new_layer.name = name;
            break;
        }
    }
    if(!is_valid_layer_name){
        std::cerr << "Error: provided layer name (" << name << ") is invalid (add_layer)" << std::endl;
    }

    //Add gas species
    new_layer.number_of_gas_species = 0;
    for(int i = 0; i < gas_species.size(); i++){
        //Iterate over all gas species: check if the name is valid before adding them to the new layer
        bool is_valid_species_name = 0;
        for(int j = 0; j < this->valid_gas_species.size(); j++){
            if(this->valid_gas_species[j] == gas_species[i]){
                is_valid_species_name = 1;
                new_layer.gas_species.push_back(gas_species[i]);
                new_layer.number_of_gas_species += 1;
                break;
            }
        }
        if(!is_valid_species_name){
            std::cerr << "Error: provided species name (" << gas_species[i] << ") is invalid (add_layer)" << std::endl;
        }
    }

    new_layer.is_dense = 0;
    new_layer.thickness = thickness;
    new_layer.nodes = nodes;
    new_layer.porosity = porosity;
    new_layer.tortuosity = tortuosity;
    new_layer.particle_diameter = particle_diameter;

    new_layer.charge_carriers = this->layer_charge_carriers.at(name);
    new_layer.number_of_charge_carriers = new_layer.charge_carriers.first + new_layer.charge_carriers.second;

    //Generate grid
    std::vector<double> center(nodes);
    std::vector<double> face(nodes+1);
    if(this->layer_count == 0){
        this->generate_grid(0.0,nodes,thickness,center,face);
        new_layer.x_center = center;
        new_layer.x_face = face;        
    }
    else{
        const int nodes_previous = this->layers[this->layer_count-1].nodes;
        const double x_begin = this->layers[this->layer_count-1].x_face[nodes_previous]; //Last face values of previous layer
        this->generate_grid(x_begin,nodes,thickness,center,face);
        new_layer.x_center = center;
        new_layer.x_face = face;        
    }


    this->layers.push_back(new_layer);
    this->layer_count++;
};


void CellConfiguration::add_gas_compartment(const bool is_air, const double inlet_pressure, const double inlet_flow_rate, 
                                            const std::vector<std::string>& gas_species, const std::vector<double>& inlet_mole_fractions){

    GasCompartment new_gas_compartment;
    new_gas_compartment.inlet_pressure = inlet_pressure;
    new_gas_compartment.inlet_flow_rate = inlet_flow_rate;

    //Check that gas_species and inlet_mole_fractions have the same length
    if(gas_species.size() != inlet_mole_fractions.size()){
        std::cerr << "Error: the length of gas_species does not equal the length of inlet_mole_fractions (is_air: " << is_air << ") (add_gas_compartment)" << std::endl;
    }
    
    //Check that inlet_mole_fractions sums to one
    double mole_fraction_sum = 0;
    for(int i = 0; i < inlet_mole_fractions.size(); i++){
        mole_fraction_sum += inlet_mole_fractions[i];
    }
    if(std::abs(mole_fraction_sum - 1.0) > 1e-5){
        std::cerr << "Error: the sum of inlet_mole_fractions does not equal one (is_air: " << is_air << ") (add_gas_compartment)" << std::endl;
    }

    //Add gas species
    new_gas_compartment.number_of_gas_species = 0;
    for(int i = 0; i < gas_species.size(); i++){
        //Iterate over all gas species: check if the name is valid before adding them to the new layer
        bool is_valid_species_name = 0;
        for(int j = 0; j < this->valid_gas_species.size(); j++){
            if(this->valid_gas_species[j] == gas_species[i]){
                is_valid_species_name = 1;
                new_gas_compartment.gas_species.push_back(gas_species[i]);
                new_gas_compartment.inlet_mole_fractions.push_back(inlet_mole_fractions[i]);
                new_gas_compartment.number_of_gas_species += 1;
                break;
            }
        }
        if(!is_valid_species_name){
            std::cerr << "Error: provided species name (" << gas_species[i] << ") is invalid (add_gas_compartment)" << std::endl;
        }
    }

    //Add new_gas_compartment to the gas_compartments pair
    if(is_air){this->gas_compartments.first = new_gas_compartment;}
    else{this->gas_compartments.second = new_gas_compartment;}
};


void CellConfiguration::check_compatibility_gas_interfaces(){


    //Identify the dense section of the MEA
    int dense_section_start_index = -1;
    int dense_section_end_index = -1;
    for(int i = 0; i < this->layer_count; i++){ //Iterate over all layers
        if(this->layers[i].is_dense){
            if(dense_section_start_index < 0){
                dense_section_start_index = i;
            }
            dense_section_end_index = i;
        }
    }

    //Check that a dense section exists
    if(dense_section_start_index < 0){
        std::cerr << "Error: no dense section of the MEA identified" << std::endl;
    }

    //Check that there are no porous layers inside the dense section
    for(int i = dense_section_start_index + 1; i < dense_section_end_index; i++){
        if(!this->layers[i].is_dense){
            std::cerr << "Error: porous layer within the dense section of the MEA" << std::endl;
        }
    }

    //Assign porous layers to the air or fuel electrode
    for(int i = 0; i < dense_section_start_index; i++){
        this->air_electrode_layer_indices.push_back(i);
        this->layers[i].is_air_electrode = 1;
    }
    for(int i = dense_section_end_index + 1; i < this->layer_count; i++){
        this->fuel_electrode_layer_indices.push_back(i);
        this->layers[i].is_air_electrode = 0;
    }

    //Check that there is at least one air electrode layer and one fuel electrode layer
    if(this->air_electrode_layer_indices.size() < 1){
        std::cerr << "Error: not at least one air electrode layer identified!" << std::endl;
    }
    if(this->fuel_electrode_layer_indices.size() < 1){
        std::cerr << "Error: not at least one fuel electrode layer identified!" << std::endl;
    }

    //Check that all layers in the air electrode have the same number of species as the air compartment
    for(int i = 0; i < this->air_electrode_layer_indices.size(); i++){
        const int layer_index = this->air_electrode_layer_indices[i];
        if(this->layers[layer_index].number_of_gas_species != this->gas_compartments.first.number_of_gas_species){
            std::cerr << "Error: porous layer (index: " << layer_index << ") within the air electrode has a different number of gas species than the air compartment" << std::endl;
        }
    }

    //Check that all porous layers in the air electrode have the same gas species as the air gas compartment, in the same order    
    for(int i = 0; i < this->gas_compartments.first.number_of_gas_species; i++){ //Iterate over all gas species in the air compartment
        for(int j = 0; j < this->air_electrode_layer_indices.size(); j++){ //Iterate over all layers in the air electrode
            const int layer_index = this->air_electrode_layer_indices[j];
            if(this->gas_compartments.first.gas_species[i] != this->layers[layer_index].gas_species[i]){
                std::cerr << "Warning: in layer " << layer_index << " (air electrode), species " << i << " is " << this->layers[layer_index].gas_species[i] << 
                ", while in the air compartment species " << i << " is " << this->gas_compartments.first.gas_species[i] << std::endl;
            }
        }
    }

    //Check that all layers in the fuel electrode have the same number of species as the fuel compartment
    for(int i = 0; i < this->fuel_electrode_layer_indices.size(); i++){
        const int layer_index = this->fuel_electrode_layer_indices[i];
        if(this->layers[layer_index].number_of_gas_species != this->gas_compartments.second.number_of_gas_species){
            std::cerr << "Error: porous layer (index: " << layer_index << ") within the fuel electrode has a different number of gas species than the fuel compartment" << std::endl;
        }
    }

    //Check that all porous layers in the fuel electrode have the same gas species as the fuel gas compartment, in the same order    
    for(int i = 0; i < this->gas_compartments.second.number_of_gas_species; i++){ //Iterate over all gas species in the fuel compartment
        for(int j = 0; j < this->fuel_electrode_layer_indices.size(); j++){ //Iterate over all layers in the fuel electrode
            const int layer_index = this->fuel_electrode_layer_indices[j];
            if(this->gas_compartments.second.gas_species[i] != this->layers[layer_index].gas_species[i]){
                std::cerr << "Warning: in layer " << layer_index << " (fuel electrode), species " << i << " is " << this->layers[layer_index].gas_species[i] << 
                ", while in the fuel compartment species " << i << " is " << this->gas_compartments.second.gas_species[i] << std::endl;
            }
        }
    }
};



void CellConfiguration::generate_grid(const double x_begin, const int nodes, const double thickness, std::vector<double>& center, std::vector<double>& face) const {
    //Uniform grid
    const double delta_x = thickness/((double) nodes);

    face[0] = x_begin;
    for(int j = 1; j < nodes + 1; j++){
        face[j] = face[j-1] + delta_x;
    }
    for(int j = 0; j < nodes; j++){
        center[j] = (face[j] + face[j+1])/2.0;
    }
        
}

//?? outdated
void CellConfiguration::print_configuration() const
{
    std::cout << "Layer name\tDense\t\tThickness [m]\tNodes\t\tGas species\tCharge carriers" << std::endl;
    for(int i = 0; i < this->layer_count; i++)
    {
        std::cout << this->layers[i].name << "\t\t" <<  this->layers[i].is_dense << "\t\t" <<  this->layers[i].thickness << "\t\t" <<  this->layers[i].nodes << "\t\t";
        //Print the gas species
        if(this->layers[i].gas_species.size() > 0){ //At least one gas species
            for(int j = 0; j < this->layers[i].gas_species.size() - 1; j++){ //Print the last species separately
                std::cout <<  this->layers[i].gas_species[j] << ",";
            }
            //Print the last species separately
            std::cout <<  this->layers[i].gas_species[this->layers[i].gas_species.size() - 1] << "\t\t"; 
        }
        else{ //No gas species
            std::cout << "\t\t"; 
        }

        //Print the charge carriers
        if(this->layers[i].charge_carriers.first + this->layers[i].charge_carriers.second > 0){ //At least one charge carrier
            std::cout << this->layers[i].charge_carriers.first << "," << this->layers[i].charge_carriers.second << std::endl;
        }
        else{ //No charge carriers
            std::cout << std::endl;
        }

    }
}


void CellConfiguration::calculate_solution_vector_unit_lengths( std::vector<int>& node_layer_indices, int& solution_vector_length, 
                                                                int& number_of_variables){
    //Set total to zero
    solution_vector_length = 0;
    number_of_variables = 0;
    this->electron_index = 0;
    this->ion_index = 0;
    int pore_interface = 0;  
    //Check if any layers have been added yet
    if(this->layer_count < 1){
        std::cerr << "Error: no layers have been set yet (calculate_solution_vector_length)" << std::endl;
    }

    for(int i = 0; i < this->layer_count; i++) //Iterate over all layers
    {
        for(int j = 0; j<this->layers[i].nodes;j++){
            node_layer_indices.push_back(i);
        }
        this->layer_start_index.push_back(solution_vector_length);
        solution_vector_length += this->layers[i].nodes;

        //Calculate the number of dependent variables
        //Mass balance
        const int n_mass = this->layers[i].number_of_gas_species;
        //Charge balance
        const int n_electron = this->layers[i].charge_carriers.first;
        const int n_ion = this->layers[i].charge_carriers.second;
        if(n_mass > this->electron_index){
            this->electron_index += n_mass-this->electron_index;
            this->ion_index += n_mass-this->electron_index;
            number_of_variables += n_mass-this->electron_index;
        }
        if(this->electron_index + n_electron > this->ion_index){
            this->ion_index += this->electron_index + n_electron - this->ion_index;
            number_of_variables += this->electron_index + n_electron - this->ion_index;
        }
        if(this->ion_index + n_ion > number_of_variables){
            number_of_variables += this->ion_index + n_ion - number_of_variables;
        }
    }
};
