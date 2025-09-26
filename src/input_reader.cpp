#include "input_reader.h"

InputReader::InputReader(){};

void InputReader::construct_class(const std::string& input_file, const std::shared_ptr<CellConfiguration>& cell_configuration){

    //Open the input file //?? check if the input file exists
    std::ifstream infile;
    infile.open(input_file); 

    //Identify line numbers where each section starts
    std::string line;
    std::vector<int> start_indices;
    int idx = 1;
    while(std::getline(infile,line))
    {
        boost::trim_left(line);
        if(line[0] == "&"[0])
        {
            start_indices.push_back(idx);
        }
        idx++;
    }
    const int last_line = idx - 1; //Line number of last line

    //Extract values from each section
    for(int i = 0; i < start_indices.size(); i++){
        //Number of entries in section (excluding the section name)
        int number_of_entries;
        if(i == start_indices.size() - 1){ //Last section
            number_of_entries = last_line - start_indices[i];
        }
        else{ //All other sections
            number_of_entries = start_indices[i+1] - start_indices[i] - 1;
        }

        //Step through every line of a section
        infile.clear(); //Clear end of file tag from infile
        infile.seekg(0); //Move back to the start

        int n = 1; //Line number
        std::string section = "";
        while(std::getline(infile,line)){
            boost::trim_left(line);
            //Skip comments
            if(line[0] == "#"[0]){
                n++;
                continue;
            }

            //Check that the section name is valid
            if(n == start_indices[i]){
                section = line.substr(1,line.length()-1);
                bool section_name_valid = 0;
                for(int j = 0; j < this->section_names.size(); j++){
                    if(section == this->section_names[j]){
                        section_name_valid = 1;
                        break;
                    }
                }
                if(!section_name_valid){
                    std::cerr << "Error: invalid section name ('" << section << "')" << std::endl;
                }
            }
            
            //Extract values from each line
            if(n > start_indices[i] && n < start_indices[i] + number_of_entries + 1){ 
                if(section == "cell_configuration"){ //For the cell_configuration section
                    //Split line into pieces
                    std::vector<std::string> pieces;
                    boost::split(pieces, line, boost::is_any_of(";"));
                    //Check that the number of pieces is correct
                    if(pieces.size() != 8)                    {
                        std::cerr << "Error: the line for " << pieces[0] << " in section &" << section << " contains " << pieces.size() << " pieces instead of 8" << std::endl;
                    }

                    //Trim pieces
                    for(int k = 1; k < pieces.size(); k++){ //Start at 1: first pieces already trimmed
                        boost::trim_left(pieces[k]);
                    }

                    //Extract information from pieces and cast to correct type
                    const std::string layer_name = pieces[0];
                    const bool is_dense = boost::lexical_cast<bool>(pieces[1]);
                    const double thickness = boost::lexical_cast<double>(pieces[2]);
                    const int nodes = boost::lexical_cast<int>(pieces[3]);

                    //Add layer
                    if(is_dense){
                        cell_configuration->add_layer(layer_name,thickness,nodes);
                    }
                    else{
                        const double porosity = boost::lexical_cast<double>(pieces[4]);
                        const double tortuosity = boost::lexical_cast<double>(pieces[5]);
                        const double particle_diameter = boost::lexical_cast<double>(pieces[6]);                        
                        std::vector<std::string> gas_species;
                        boost::split(gas_species, pieces[7], boost::is_any_of(",")); //Split the last piece into the separate gas species                        
                        cell_configuration->add_layer(layer_name,thickness,nodes,porosity,tortuosity,particle_diameter,gas_species);
                    }
                }
                else if(section == "gas_compartment"){ //For the gas compartment section
                    //Split line into pieces
                    std::vector<std::string> pieces;
                    boost::split(pieces, line, boost::is_any_of(";"));
                    //Check that the number of pieces is correct
                    if(pieces.size() != 5){
                        std::cerr << "Error: the line for " << pieces[0] << " (in section &" << section << ") contains " << pieces.size() << " pieces instead of 5" << std::endl;
                    }

                    //Trim pieces
                    for(int k = 1; k < pieces.size(); k++){ //Start at 1: first pieces already trimmed
                        boost::trim_left(pieces[k]);
                    }

                    const std::string compartment_name = pieces[0];
                    const bool is_air = compartment_name == "air" ? 1 : 0;
                    const double inlet_pressure = boost::lexical_cast<double>(pieces[2]);
                    const double inlet_flow_rate = boost::lexical_cast<double>(pieces[4]);
                    std::vector<std::string> gas_species;
                    boost::split(gas_species, pieces[1], boost::is_any_of(","));
                    std::vector<std::string> inlet_mole_fractions_string;
                    boost::split(inlet_mole_fractions_string, pieces[3], boost::is_any_of(","));
                    std::vector<double> inlet_mole_fractions;
                    for(int k = 0; k < inlet_mole_fractions_string.size(); k++){
                        inlet_mole_fractions.push_back(boost::lexical_cast<double>(inlet_mole_fractions_string[k]));
                    }

                    //Add gas compartment
                    cell_configuration->add_gas_compartment(is_air,inlet_pressure,inlet_flow_rate,gas_species,inlet_mole_fractions);
                }
                else{   //For all remaining sections
                    //Split line into pieces
                    std::vector<std::string> pieces;
                    boost::split(pieces, line, boost::is_any_of(";"));
                    //Check that the number of pieces is correct
                    if(pieces.size() != 4){
                        std::cerr << "Error: the line for " << pieces[0] << " contains " << pieces.size() << " pieces instead of 4" << std::endl;
                    }

                    //Trim pieces
                    for(int k = 1; k < pieces.size(); k++){ //Start at 1: first pieces already trimmed
                        boost::trim_left(pieces[k]);
                    }

                    boost::variant<int,bool,double,std::string> value = 0;
                    //Check for data type, cast to correct data type and assign to boost::variant type variable for storage
                    if(pieces[2] == "string"){
                        value = boost::lexical_cast<std::string>(pieces[1]);
                    }
                    else if(pieces[2] == "bool"){
                        value = boost::lexical_cast<bool>(pieces[1]);
                    }
                    else if(pieces[2] == "int"){
                        value = boost::lexical_cast<int>(pieces[1]);
                    }
                    else if(pieces[2] == "double"){
                        value = boost::lexical_cast<double>(pieces[1]);
                    }
                    else{
                        std::cerr << "Error: invalid data type for " << pieces[0] << ", use string, bool, int or double" << std::endl;
                    }
                    //Store pieces
                    this->input_values.insert(std::make_pair(pieces[0],std::make_pair(value,pieces[2]))); //key + (value + type)
                }

            }
        
            n++;
        }
    }
    infile.close();

    cell_configuration->check_compatibility_gas_interfaces();

}


const int InputReader::get_int_input(const std::string& key) const
{
    std::pair<boost::variant<int,bool,double,std::string>,std::string> output;

    if(this->input_values.count(key) == 1) //Check if key is present in input_values
    {
        output = this->input_values.at(key);
        if(output.second != "int")
        {
            std::cerr << "Error: wrong get method chosen for " << key << ": get_int instead of get_" << output.second << std::endl;
        }   
    }
    else{
        std::cerr << "Error: value requested for " << key << " not present in unordered map (get int)" << std::endl;
    }
    return boost::get<int>(output.first);    
};

const bool InputReader::get_bool_input(const std::string& key) const
{
    std::pair<boost::variant<int,bool,double,std::string>,std::string> output;

    if(this->input_values.count(key) == 1) //Check if key is present in input_values
    {
        output = this->input_values.at(key);
        if(output.second != "bool")
        {
            std::cerr << "Error: wrong get method chosen for " << key << ": get_bool instead of get_" << output.second << std::endl;
        }   
    }
    else{
        std::cerr << "Error: value requested for " << key << " not present in unordered map (get bool)" << std::endl;
    }
    return boost::get<bool>(output.first);    
};

const double InputReader::get_double_input(const std::string& key) const
{
    std::pair<boost::variant<int,bool,double,std::string>,std::string> output;

    if(this->input_values.count(key) == 1) //Check if key is present in input_values
    {
        output = this->input_values.at(key);
        if(output.second != "double")
        {
            std::cerr << "Error: wrong get method chosen for " << key << ": get_double instead of get_" << output.second << std::endl;
        }   
    }
    else{
        std::cerr << "Error: value requested for " << key << " not present in unordered map (get double)" << std::endl;
    }
    return boost::get<double>(output.first);    
};

const std::string InputReader::get_string_input(const std::string& key) const
{
    std::pair<boost::variant<int,bool,double,std::string>,std::string> output;

    if(this->input_values.count(key) == 1) //Check if key is present in input_values
    {
        output = this->input_values.at(key);
        if(output.second != "string")
        {
            std::cerr << "Error: wrong get method chosen for " << key << ": get_string instead of get_" << output.second << std::endl;
        }   
    }
    else{
        std::cerr << "Error: value requested for " << key << " not present in unordered map (get string)" << std::endl;
    }
    return boost::get<std::string>(output.first);    
};