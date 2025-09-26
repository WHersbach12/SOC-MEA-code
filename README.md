# SOC-MEA-code
# Purpose 
This code resolves the profiles for gasses and charge carriers inside a membrane electrode assembly in 1D for a solid oxide cell. It was made by Wisse Hersbach under the supervision of Aayan Banerjee. The code allows for multiple different layers to be configured for the solid oxide cell with different microstructures, which allows for fast estimation of the performance of a single cell.

# Installation
This installation guide is for Debian 12 (or Windows, using the Debian WSL from the Microsoft Store)

## Dependencies

### Step I: Debian version
To ensure the correct version of PETSc is installed in a later step, Debian 12 should be used. Check this using the following command:
```
cat /etc/issue
```

### Step II: Update and install essential dependencies
```
sudo apt-get update
sudo apt-get install -y build-essential git-all cmake libopenmpi-dev libboost-all-dev pkg-config wget gfortran python3
```
### Step III: Update and install python packages
```
sudo apt-get install python3-matplotlib
```
### Step IV: Install PETSc
```
sudo apt install petsc-dev
```
Note that the makers of PETSc also provide an alternative installation method for PETSc. However, that method was found to be significantly harder to combine correctly with PkgConfig in CMake. As such, the above is recommended - although this method limits the available PETSc versions to the PETSc source package versions provided by Debian.

# IntelliSense
If the Visual Studio Code IDE is used together with the Debian WSL, add "/usr/include/petsc" to includePath in c_cpp_properties.json and set compilerPath to "/usr/bin/mpiCC" for the inline compiler to work correctly.

# Building 
Enter the directory and create a build directory:
```
cd button-cell
mkdir build
cd build
```
Call CMake:
```
cmake ..
```
Compile:
```
make -jN
```
where N is the number of threads to be used in parallel for the compilation.

# Execution
```
mpirun -n N ./nernst -i ../input/ -o ../output/
```
where N is the number of cores with which to execute the program. Note that the user is free to change the relative input and output paths to different directories. All files in the chosen output directory will be deleted and replaced by the output of this program.
