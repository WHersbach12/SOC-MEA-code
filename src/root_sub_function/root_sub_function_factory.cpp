#include "root_sub_function_factory.h"


RootSubFunctionFactory::RootSubFunctionFactory(){

};


std::shared_ptr<RootSubFunction> RootSubFunctionFactory::factory_method(const int rank, const double temperature, const std::string& layer_name, 
                                                                        const std::vector<double>& x_center, const std::vector<double>& x_face,
                                                                        const std::vector<std::string>& gas_species, const bool galvanostatic,
                                                                        const double applied_I_or_V, const double area, const int layer_index, 
                                                                        const std::shared_ptr<CellConfiguration>& cell_configuration) const{
    
    std::shared_ptr<RootSubFunction> output;

    if(layer_name == "LSCF"){
        output = std::make_shared<RootSubFunctionLSCF>( rank,x_center,x_face,gas_species,temperature,area,galvanostatic,
                                                        applied_I_or_V,layer_index,cell_configuration);
    }
    else if(layer_name == "Ni-YSZ"){
        output = std::make_shared<RootSubFunctionNiYSZ>(rank,x_center,x_face,gas_species,temperature, galvanostatic,
                                                        applied_I_or_V,area,layer_index,cell_configuration);
    }    
    else if(layer_name == "Ni-CGO"){
        output = std::make_shared<RootSubFunctionNiCGO>(rank,x_center,x_face,gas_species,temperature,area,galvanostatic,
                                                        applied_I_or_V,layer_index,cell_configuration);
    }        
    else if(layer_name == "YSZ"){
        output = std::make_shared<RootSubFunctionYSZ>(rank,x_center,x_face,temperature,layer_index,cell_configuration);
    }    
    else if(layer_name == "3YSZ"){
        output = std::make_shared<RootSubFunction3YSZ>(rank,x_center,x_face,temperature,layer_index,cell_configuration);
    }   
    else if(layer_name == "CeSSZ"){
        output = std::make_shared<RootSubFunctionCeSSZ>(rank,x_center,x_face,temperature,layer_index,cell_configuration);
    }   
    else if(layer_name == "Support"){
        output = std::make_shared<RootSubFunctionMetalSupport>(rank,x_center,x_face,gas_species,temperature,area,galvanostatic,
                                                        applied_I_or_V,layer_index,cell_configuration);
    }   
    else if(layer_name == "CGO"){
        output = std::make_shared<RootSubFunctionCGO>(rank,x_center,x_face,gas_species,temperature,area,galvanostatic,
                                                        applied_I_or_V,layer_index,cell_configuration);
    }    
    else if(layer_name == "LSM-YSZ"){
        output = std::make_shared<RootSubFunctionLSMYSZ>(rank,x_center,x_face,gas_species,temperature,area,galvanostatic,
                                                        applied_I_or_V,layer_index,cell_configuration);
    }           
    else if(layer_name == "LSM"){
        output = std::make_shared<RootSubFunctionLSM>(rank,x_center,x_face,gas_species,temperature,area,galvanostatic,
                                                        applied_I_or_V,layer_index,cell_configuration);
    }           
    else{
        std::cerr << "Error: provided layer name is invalid (RootSubFunctionFactory::factory_method)" << std::endl;
    }

    return output;
};