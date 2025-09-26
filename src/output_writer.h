#include "cell_configuration.h"

class OutputWriter{
private:
    const int precision = 8; //Floating point writing precision

    int rank;
    std::string output_folder;

public: 
    OutputWriter(const int _rank, const std::string& _output_folder, const std::string& input_folder);

    // PetscErrorCode write_solution_to_csv(   const int* nodes_per_process, const int comm_size, const int solution_vector_length, 
    //                                         const std::shared_ptr<CellConfiguration>& cell_configuration, Vec x) const;
    //Creates a csv file with two columns: x-coordinate - solution vector
    PetscErrorCode write_solution_to_csv(   const int* nodes_per_process, const int comm_size, const int solution_vector_length, 
                                            const std::shared_ptr<CellConfiguration>& cell_configuration, Vec x, Vec Fluxes, Vec Source_term, DM da) const;

private:




};