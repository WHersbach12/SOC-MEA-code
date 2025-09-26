import matplotlib.pyplot as plt
import numpy as np
import math
import os
import csv


#?? update input reader

F = 96485.33 #Faraday's constant
R = 8.31446
mu_H2O_0 = -237140; #https://www.emse.fr/~bguy/textes%20pdf%20etc/albert%20guy%20damidot.pdf
rgb = np.array([[33,183,130], [101, 41, 231], [41,140,231], [28,39,50], [20,177,188], [46, 64, 217], [162, 23, 213]])
rgb = rgb/250
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=rgb)

#Define the Layer and Layers class
class Layer:
    def __init__(self,name,is_dense,thickness,nodes,gas_species,charge_carriers):
        self.name = name
        self.is_dense = is_dense
        self.thickness = thickness
        self.nodes = nodes
        self.gas_species = gas_species
        self.charge_carriers = charge_carriers

        #Count the number of gas species
        self.number_of_gas_species = 0
        if len(self.gas_species) == 1:
            if self.gas_species[0] == "-":
                self.number_of_gas_species = 0
            else:
                self.number_of_gas_species = 1
        else:
            self.number_of_gas_species = len(self.gas_species)

        #Count the number of charge carriers
        self.number_of_charge_carriers = self.charge_carriers[0] + self.charge_carriers[1]

        #Calculate the number of dependent variables
        self.variables = self.number_of_gas_species + self.number_of_charge_carriers


class Layers:
    def __init__(self):
        self.layers = []
        self.layer_count = 0

        #Create the dictionary describing the charge carriers present in each type of layer
        #(also defines the valid layer types)
        self.layer_charge_carriers = {  "LSCF":(1,1),
                                        "LSCF-CGO":(1,2),
                                        "LSM-YSZ":(1,1),
                                        "YSZ":(0,1),
                                        "3YSZ":(0,1),
                                        "Ni-YSZ":(1,1),
                                        "Ni-CGO":(1,1),
                                        "Support":(1,0),
                                        "CGO":(1,1),
                                        "LSM":(1,0)}



    def add_layer(self,name,is_dense,thickness,nodes,gas_species):
        charge_carriers = self.layer_charge_carriers[name]
        layer = Layer(name,is_dense,thickness,nodes,gas_species,charge_carriers)
        self.layers.append(layer)
        self.layer_count += 1





#?? read input arguments to call

output_folder = "../output/"


#Read input file
input_values = {} #Dictionary to contain the input values
cell_configuration = Layers() #Class to contain the cell configuration

input_file = output_folder + "input/" + "input.txt"
data = np.genfromtxt(input_file, dtype = str, delimiter="\n")
rows = data.shape[0]

#Identify line numbers where each section starts
section_start_indices = []
for i in range(0,rows):
    if data[i][0] == "&":
        section_start_indices.append(i)

#Extract values from each section
for i in range(0,len(section_start_indices)): #Iterate over all sections

    #Identify the starting index of the next section
    section_end = -1
    if i < len(section_start_indices) - 1:
        section_end = section_start_indices[i+1]
    else:
        section_end = rows
    
    if data[section_start_indices[i]][1:] == "cell_configuration": #Cell configuration section
        for j in range(section_start_indices[i]+1,section_end): #Iterate over all rows of this section
            pieces = data[j].split(";") 
            for k in range(0,len(pieces)):
                pieces[k] = pieces[k].strip()

            #Create list of gas species
            gas_species = pieces[-1].split(",")

            #Add layer
            cell_configuration.add_layer(pieces[0],bool(pieces[1]=='1'),float(pieces[2]),int(pieces[3]),gas_species)


    else: #All other sections
        for j in range(section_start_indices[i]+1,section_end): #Iterate over all rows of this section
            pieces = data[j].split(";") 
            for k in range(0,len(pieces)):
                pieces[k] = pieces[k].strip()
            #Add entry to input_values
            input_values[pieces[0]] = pieces[1]

T = float(input_values["T"])
dH_O2 = 0 + 31.32234*T/1000 - 20.23531*pow(T/1000,2)/2 + 57.8664*pow(T/1000,3)/3 - 36.50624*pow(T/1000,4)/4 - 0.007374*pow(T/1000,-1)
dS_O2 = 205.152
mu_O2_0 = dH_O2 - T*dS_O2
dH_H2 = 0 + 33.066178*T/1000 - 11.3634*pow(T/1000,2)/2 + 11.4328*pow(T/1000,3)/3 - 2.772874*pow(T/1000,4)/4 - 0.158558*pow(T/1000,-1)
dS_H2 = 130.68
dH_H2O = -241.83e3 + 30.092*T/1000 - 6.832514*pow(T/1000,2)/2 + 6.79345*pow(T/1000,3)/3 - 2.53448*pow(T/1000,4)/4 + 0.082139*pow(T/1000,-1)
dS_H2O = 188.853
mu_H2_0 = dH_H2 - T*dS_H2
mu_H2O_0 = dH_H2O - T*dS_H2O

#Extract the solution vector from file
solution_vector = np.loadtxt(output_folder + "solution_vector.txt", dtype = float,delimiter = ' ')
#Determine the start and end index for every layer
layer_start_index = np.zeros(cell_configuration.layer_count,dtype=np.int64) #Start index for every layer
layer_end_index = np.zeros(cell_configuration.layer_count,dtype=np.int64) #End index for every layer
for i in range(1,cell_configuration.layer_count): #Start from the second layer
    layer_start_index[i] = layer_start_index[i-1] + cell_configuration.layers[i-1].nodes
    layer_end_index[i-1] = layer_start_index[i]
layer_end_index[-1] = len(solution_vector)

Air_electrode = 1
Fuel_electrode = 0
electron_index = cell_configuration.layers[-1].number_of_gas_species + 1
ion_index = cell_configuration.layers[-1].number_of_gas_species + cell_configuration.layers[-1].number_of_charge_carriers

p_air =  np.empty([1, cell_configuration.layers[0].number_of_gas_species+1])
p_fuel = np.empty([1, cell_configuration.layers[-1].number_of_gas_species+1])
j_air =  np.empty([1, cell_configuration.layers[0].number_of_gas_species+1])
j_fuel = np.empty([1, cell_configuration.layers[-1].number_of_gas_species+1])
mu_e = solution_vector[:,[0, electron_index]]
mu_o = solution_vector[:,[0, ion_index]]
source_e = solution_vector[:,[0, ion_index + electron_index]] 
j_e = solution_vector[:,[2*ion_index+1, 2*ion_index+1 + electron_index]]
j_o = solution_vector[:,[2*ion_index+1, 2*ion_index+1 + ion_index]]
chi = -(-2*mu_e[:,1] + mu_o[:,1] + 0.5*mu_O2_0)/(2.0*F)
eta = (2*mu_e[:,1] - mu_o[:,1] + (mu_H2_0 - mu_H2O_0))/(2.0*F)
c = 1e5/R/T
#Separate the solution vector in the important bits
for n in range(cell_configuration.layer_count):
    if cell_configuration.layers[n].is_dense:
        if cell_configuration.layers[n].charge_carriers[0] == 0: #If layer doesn't conduct electrons, dont plot the values
            mu_e[layer_start_index[n]:layer_end_index[n],1] = np.nan
            source_e[layer_start_index[n]:layer_end_index[n],1] = np.nan
            eta[layer_start_index[n]:layer_end_index[n]] = np.nan
            chi[layer_start_index[n]:layer_end_index[n]] = np.nan
            j_e[layer_start_index[n]:layer_end_index[n],1] = np.nan
        Air_electrode = 0
        Fuel_electrode = 1

    elif Air_electrode:
        p_air = np.append(p_air,solution_vector[layer_start_index[n]:layer_end_index[n], 0:(cell_configuration.layers[n].number_of_gas_species+1)],axis=0) 
        j_air = np.append(j_air,solution_vector[layer_start_index[n]:layer_end_index[n]+1, 2*ion_index+1:(2*ion_index+1 + cell_configuration.layers[n].number_of_gas_species+1)],axis=0) 
        eta[layer_start_index[n]:layer_end_index[n]] = np.nan
        chi[layer_start_index[n]:layer_end_index[n]] -= 0.5*R*T*np.log(solution_vector[layer_start_index[n]:layer_end_index[n],1]/(solution_vector[layer_start_index[n]:layer_end_index[n],1]+solution_vector[layer_start_index[n]:layer_end_index[n],2]))/(2.0*F)
        if cell_configuration.layers[n].charge_carriers[0] == 0: #If layer doesn't conduct electrons, dont plot the values
            mu_e[layer_start_index[n]:layer_end_index[n],1] = np.nan
            chi[layer_start_index[n]:layer_end_index[n]] = 0.0
        if cell_configuration.layers[n].charge_carriers[1] == 0: #If layer doesn't conduct oxygen anions, dont plot the values
            mu_o[layer_start_index[n]:layer_end_index[n],1] = np.nan
            chi[layer_start_index[n]:layer_end_index[n]] = 0.0

    elif Fuel_electrode:
        p_fuel = np.append(p_fuel,solution_vector[layer_start_index[n]:layer_end_index[n], 0:(cell_configuration.layers[n].number_of_gas_species+1)],axis=0) 
        j_fuel = np.append(j_fuel,solution_vector[layer_start_index[n]:layer_end_index[n], 2*ion_index+1:(2*ion_index+1 + cell_configuration.layers[n].number_of_gas_species+1)],axis=0) 
        chi[layer_start_index[n]:layer_end_index[n]] = np.nan
        eta[layer_start_index[n]:layer_end_index[n]] +=  R*T*np.log(solution_vector[layer_start_index[n]:layer_end_index[n],1]/solution_vector[layer_start_index[n]:layer_end_index[n],2])/(2.0*F)
        if cell_configuration.layers[n].charge_carriers[0] == 0: #If layer doesn't conduct electrons, dont plot the values
            mu_e[layer_start_index[n]:layer_end_index[n],1] = np.nan
            j_e[layer_start_index[n]:layer_end_index[n],1] = np.nan

        if cell_configuration.layers[n].charge_carriers[1] == 0: #If layer doesn't conduct oxygen anions, dont plot the values
            mu_o[layer_start_index[n]:layer_end_index[n],1] = np.nan
            j_o[layer_start_index[n]:layer_end_index[n],1] = np.nan

#Create figure objects
fig_mass, ax_mass               = plt.subplots(figsize=(6,4.5))
fig_charge, ax_charge           = plt.subplots(figsize=(6,4.5))
fig_mass_flux, ax_mass_flux     = plt.subplots(figsize=(6,4.5))
fig_charge_flux, ax_charge_flux = plt.subplots(figsize=(6,4.5))
fig_source, ax_source           = plt.subplots(figsize=(6,4.5))
fig_overPot, ax_overPot         = plt.subplots(figsize=(6,4.5))

#Setup charge balance figure
#Add vertical lines between layers
layer_right_boundary = 0
for i in range(0,cell_configuration.layer_count-1): 
    layer_right_boundary += cell_configuration.layers[i].thickness*1e6
    ax_charge.axvline(x=layer_right_boundary, color = 'k', linewidth = 1) #For charge balance
    ax_mass.axvline(x=layer_right_boundary, color = 'k', linewidth = 1) #For mass balance
    ax_mass_flux.axvline(x=layer_right_boundary, color = 'k', linewidth = 1) #For mass balance
    ax_charge_flux.axvline(x=layer_right_boundary, color = 'k', linewidth = 1) #For mass balance
    ax_overPot.axvline(x=layer_right_boundary, color = 'k', linewidth = 1) #For mass balance
    ax_source.axvline(x=layer_right_boundary, color = 'k', linewidth = 1) #For mass balance
#Add text for layer material
for i in range(0,cell_configuration.layer_count):
    layer_center_index = layer_start_index[i] + round((layer_end_index[i]-layer_start_index[i])/2)-1
    layer_center_position = solution_vector[layer_center_index,0] 
    layer_center_frac = layer_center_position/solution_vector[-1,0] - 0.03*(layer_end_index[i]-layer_start_index[i])/layer_center_index
    yloc = 0.95
    ax_charge.annotate(cell_configuration.layers[i].name,xy=(layer_center_frac,yloc),xycoords='axes fraction', xytext=(layer_center_frac,yloc),textcoords='axes fraction')
    ax_mass.annotate(cell_configuration.layers[i].name,xy=(layer_center_frac,yloc),xycoords='axes fraction', xytext=(layer_center_frac,yloc),textcoords='axes fraction')
    ax_source.annotate(cell_configuration.layers[i].name,xy=(layer_center_frac,yloc),xycoords='axes fraction', xytext=(layer_center_frac,yloc),textcoords='axes fraction')
    ax_mass_flux.annotate(cell_configuration.layers[i].name,xy=(layer_center_frac,yloc),xycoords='axes fraction', xytext=(layer_center_frac,yloc),textcoords='axes fraction')
    ax_charge_flux.annotate(cell_configuration.layers[i].name,xy=(layer_center_frac,yloc),xycoords='axes fraction', xytext=(layer_center_frac,yloc),textcoords='axes fraction')
    ax_overPot.annotate(cell_configuration.layers[i].name,xy=(layer_center_frac,yloc),xycoords='axes fraction', xytext=(layer_center_frac,yloc),textcoords='axes fraction')

#Setup charge balance figure
plt.figure(fig_charge)
ax_charge.plot(mu_e[:,0]*1e6,mu_e[:,1]/F,label=r"$\tilde{\mu}_{\mathrm{e^-}}^{*}$")
ax_charge.plot(mu_o[:,0]*1e6,mu_o[:,1]/(2.0*F),label=r"$\tilde{\mu}_{\mathrm{O^{2-}}}^{*}$")
ax_charge.legend()
ax_charge.set_xlabel("$x$ $\mathrm{[\mu m]}$")
ax_charge.set_ylabel(r"$\tilde{\mu}^{*}$ $\mathrm{[V]}$")
plt.savefig(output_folder + "Charge_balance" +  ".PNG", bbox_inches = 'tight')
plt.close()

#Setup mass balance figure
plt.figure(fig_mass)
for j in range(0,cell_configuration.layers[0].number_of_gas_species):    
    #plot data
    mass_label = 'Air side '+cell_configuration.layers[0].gas_species[j]
    ax_mass.plot(p_air[1:,0]*1e6, p_air[1:,j+1]/(np.sum(p_air[1:,1:],1)), label = mass_label)
    ax_mass_flux.plot(j_air[1:,0]*1e6, j_air[1:,j+1], label = mass_label)

for j in range(0,cell_configuration.layers[-1].number_of_gas_species):    
    #plot data
    mass_label = 'Fuel side '+cell_configuration.layers[-1].gas_species[j]
    ax_mass.plot(p_fuel[1:,0]*1e6,p_fuel[1:,j+1]/(np.sum(p_fuel[1:,1:],1)), label = mass_label)
    ax_mass_flux.plot(j_fuel[1:,0]*1e6,j_fuel[1:,j+1], label = mass_label)
    
ax_mass.set_xlabel("$x$ $\mathrm{[\mu m]}$")
ax_mass.set_ylabel("$x_i$ $[-]$")
ax_mass.legend()
plt.savefig(output_folder + "Mass_balance" + ".PNG", bbox_inches = 'tight')
plt.close()

#Setup mass fluxes figure
plt.figure(fig_mass_flux)
ax_mass_flux.set_xlabel("$x$ $\mathrm{[\mu m]}$")
ax_mass_flux.set_ylabel("$J_i$ $[mol/m^2/s]$")
ax_mass_flux.legend()
plt.savefig(output_folder + "Mass_flux" + ".PNG", bbox_inches = 'tight')
plt.close()

#Setup charge fluxes figure
plt.figure(fig_charge_flux)
ax_charge_flux.plot(j_e[:,0]*1e6,j_e[:,1]/1e4,label=r"$J_{\mathrm{e^-}}$")
ax_charge_flux.plot(j_o[:,0]*1e6,j_o[:,1]/1e4,label=r"$J_{\mathrm{O^{2-}}}$")
ax_charge_flux.set_xlabel("$x$ $\mathrm{[\mu m]}$")
ax_charge_flux.set_ylabel("$J_i$ $[A/cm^2]$")
ax_charge_flux.legend()
plt.savefig(output_folder + "Charge_flux" + ".PNG", bbox_inches = 'tight')
plt.close()

#Setup sources figure
plt.figure(fig_source)
ax_source.plot(source_e[:,0]*1e6,source_e[:,1], label = r"$i^v_e$")
ax_source.set_xlabel("$x$ $\mathrm{[\mu m]}$")
ax_source.set_ylabel("$i$ $[A/m^3]$")
ax_source.legend()
plt.savefig(output_folder + "Current_density" + ".PNG", bbox_inches = 'tight')
plt.close()

#Setup over potential figure
plt.figure(fig_overPot)
ax_overPot.plot(mu_e[:,0]*1e6,eta, label=r"$\eta^fe$")
ax_overPot.plot(mu_e[:,0]*1e6,chi, label=r"$\eta^ae$")
ax_overPot.set_xlabel("$x$ $\mathrm{[\mu m]}$")
ax_overPot.set_ylabel("$\eta$ $[V]$")
ax_overPot.legend()
plt.savefig(output_folder + "Overpotential" + ".PNG", bbox_inches = 'tight')
plt.close()
