//Standard cpp libraries
#include <string>
#include <vector>
#include <fstream>          //Reading/writing from files
#include <memory>           //Smart pointers
#include <iostream>         //Print to terminal
#include <chrono>           //Track elapsed time
#include <unordered_map>    //Dictionaries
#include <math.h>           //For pi
#include <filesystem>       //For checking for the existence/creation of files, folders

//General boost libraries
#include <boost/mpi.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/variant.hpp> //Serialization of std::variant not possible

//Boost serialization libraries for broadcasting of complex data types
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/array.hpp>

//Petsc libraries
#include <petsc.h>
#include <petscsystypes.h>
#include <petscsnes.h>
#include <petscsys.h>


//List of used concepts, in order of decreasing importance
//. MPI (+ serialization)               Hard (https://mpitutorial.com/tutorials/mpi-send-and-receive/)
//. const keyword (all uses)            Intermediate - not before classes    
//. Classes + encapsulation             Intermediate
//. Function overloading                Easy
//. Class inheritance and polymorphism  Intermediate - not before classes
//. Factory method design pattern       Intermediate - not before class polymorphism
//. Arrays and pointers                 Intermediate
//. Ternary operator                    Easy