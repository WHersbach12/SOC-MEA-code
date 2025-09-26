#include "root_sub_function_LSCF.h"
#include "root_sub_function_YSZ.h"
#include "root_sub_function_3YSZ.h"
#include "root_sub_function_CeSSZ.h"
#include "root_sub_function_NiCGO.h"
#include "root_sub_function_NiYSZ.h"
#include "root_sub_function_metal_support.h"
#include "root_sub_function_CGO.h"
#include "root_sub_function_LSMYSZ.h"
#include "root_sub_function_LSM.h"

class RootSubFunctionFactory{
private:


public:
    RootSubFunctionFactory();

    std::shared_ptr<RootSubFunction> factory_method(const int rank, const double temperature, const std::string& layer_name, 
                                                    const std::vector<double>& x_center, const std::vector<double>& x_face,
                                                    const std::vector<std::string>& gas_species, const bool galvanostatic,
                                                    const double applied_I_or_V, const double area,
                                                    const int layer_index, const std::shared_ptr<CellConfiguration>& cell_configuration) const;

private:


};

