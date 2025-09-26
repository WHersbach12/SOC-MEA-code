cmake_minimum_required(VERSION 3.18)  #Set minimum cmake requirements
set(CMAKE_CXX_COMPILER "/usr/bin/mpiCC") #Set compiler: MPI wrapper of g++ compiler - required for PETSc
set(PROJECT_NAME nernst)
project (${PROJECT_NAME}) #Name project

#Packages
find_package(PkgConfig) #PkgConfig 
if (${PKG_CONFIG_FOUND}) 
  pkg_check_modules(PETSC PETSc) #PETSc
  find_package(Boost 1.74 COMPONENTS mpi REQUIRED) #Boost
else()
  message(WARNING "PkgConfig not found")
endif()

#PETSc
if (${PETSC_FOUND})
  include_directories(${PETSC_INCLUDE_DIRS})
else()
  message(WARNING "PETSc not found")
endif()

#PETSc
if (${BOOST_FOUND})
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message(WARNING "Boost not found")
endif()

#Create executable
file(GLOB_RECURSE SOURCES "src/*.cpp") #Locate all source files in the /src/ folder (and it directories)
#list(REMOVE_ITEM SOURCES ..) #Command to remove specific specific file from the sources list
add_executable(${PROJECT_NAME} ${SOURCES}) #Create executable

#Set C++20
add_definitions(-std=c++20 -march=native)
add_definitions(-O3) #Compiler optimization setting

#Link libraries
target_link_libraries(${PROJECT_NAME} ${PETSC_LIBRARIES} 
                                      ${Boost_LIBRARIES})


