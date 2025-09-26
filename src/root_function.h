#include "root_sub_function/root_sub_function_factory.h"

//?? move this 
//Species index (m: order as provided in input file,
//c: electrons before oxygen anions - main conductor first)


typedef std::unordered_map<int,std::string> umpas_is;

class RootFunction
{
private:
    //Process specific properties
    int rank;

    //Vectors for the Jacobian (not broadcasted)
    Vec f;
    Vec df;

    std::shared_ptr<CellConfiguration> cell_configuration; //Cannot be broadcasted: address would change

    std::shared_ptr<RootSubFunction>* root_sub_functions; //Not broadcasted

    int* layer_index_mapping; //Not broadcasted

    umpas_is petsc_convergence_reasons;

    //Broadcasted fields
    double temperature;

    double area;

    bool galvanostatic;

    double applied_I_or_V;

    int comm_size;

    int solution_vector_length;

    int* nodes_per_process;

    std::vector<int> node_layer_index; //The layer index for every unit of the solution vector

    typedef struct{
        DM da;
        std::unique_ptr<RootFunction> rf_wrap;
    }   AppCtx;

    //Configure serialization
    friend class boost::serialization::access;

    template<class Archive> void serialize(Archive& ar, const unsigned int version){
        ar & temperature;
        ar & area;
        ar & galvanostatic;
        ar & applied_I_or_V;
        ar & comm_size;
        ar & solution_vector_length;
        ar & nodes_per_process;
        ar & node_layer_index;
    }

public:
    Vec Fluxes;
    Vec Source_term;

    //Construction and destruction
    RootFunction(const std::shared_ptr<CellConfiguration>& _cell_configuration); //Constructor

    ~RootFunction(); //Destructor

    void construct_class(   const double _temperature,
                            const double _area,
                            const bool _galvanostatic,
                            const double _applied_I_or_V,
                            const int _comm_size,
                            const int _solution_vector_length,
                            int* _nodes_per_process,
                            const std::vector<int> _node_layer_index);


    PetscErrorCode construct_class_per_process(const int _rank, Vec x, DM da); 
    //x is passed to create the necessary vectors for the Jacobian by duplication

    //Public member functions
    PetscErrorCode set_initial_guess(Vec x, DM da) const; //?? Location?

    PetscErrorCode calculate_root_function(Vec x, Vec f, DM da) const;

    PetscErrorCode calculate_jacobian(Vec x, Mat J, DM da) const;

    const std::string get_petsc_convergence_reasons(const int key) const;

    PetscErrorCode writer(Vec x, PetscScalar* xx, PetscScalar* flux, PetscScalar* source, DM da) const;

private:

 

};


PetscErrorCode root_function_wrapper(SNES, Vec, Vec, void *);

PetscErrorCode jacobian_wrapper(SNES, Vec, Mat, Mat, void *);

PetscErrorCode write_wrapper(const int, const std::string&, const std::string&, Vec x, void *);