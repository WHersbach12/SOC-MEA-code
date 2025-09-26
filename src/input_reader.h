#include "cell_configuration.h"

#pragma once

typedef std::unordered_map<std::string,std::pair<boost::variant<int,bool,double,std::string>,std::string>> umap_sv; //key + (value + type)

//Stateful
class InputReader
{
private:

    umap_sv input_values;

    std::vector<std::string> section_names = {
        "run",
        "operation",
        "cell_configuration",
        "gas_compartment"
    };
    //?? should be constant, not possible with serialization

    //Configure serialization
    friend class boost::serialization::access;

    template<class Archive> void serialize(Archive& ar, const unsigned int version){
        ar & input_values;
        ar & section_names;
    };

    // template<class Archive> void save_construct_data(Archive& ar, const InputReader* t, const unsigned int version)
    // {
    //     ar << t->section_names;
    // }

    // template<class Archive> void load_construct_data(Archive& ar, InputReader* t, const unsigned int version)
    // {
    //     new (t) InputReader();
    // }


public:
    InputReader();

    void construct_class(const std::string& input_file, const std::shared_ptr<CellConfiguration>& cell_configuration);

    //Retrieve unordered map
    inline const umap_sv& get_input_values() const{
        return this->input_values;
    }     

    const int get_int_input(const std::string& key) const;

    const bool get_bool_input(const std::string& key) const;

    const double get_double_input(const std::string& key) const;

    const std::string get_string_input(const std::string& key) const;


private:


};