#!/usr/bin/env python
# coding: utf-8

# ## Cardiac Regulatory Network for CMRG
# ### Originally developed by Shulin Cao
# ### Modified/annotated by Mihir Samdarshi
# ### Only for use in CMRG

###------Network Simulator------###
###---------Shulin Cao----------###
###-------Mihir Samdarshi-------###
###------CMRG, UC San Diego-----###

###import packages###
import pandas as pd
import collections
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import ntpath
import os
sns.set()

######################################
# SET THE EXCEL SHEET HERE #
######################################
macrophage = "/Users/mihir/Documents/Summer/Models/experimental_models/macrophage_model.xlsx"
combined_with_cpd43 = "/Users/mihir/Documents/Summer/Models/experimental_models/combined_model_with_cpd43.xlsx"
combined_without_cpd43 = "/Users/mihir/Documents/Summer/Models/experimental_models/combined_model_without_cpd43.xlsx"
og_fibroblast = "/Users/mihir/Documents/Summer/Models/original_models/original_fibroblast_model.xlsx"
og_cardiomyocyte = "/Users/mihir/Documents/Summer/Models/original_models/original_cardiomyocyte_model.xlsx"
ficks = "/Users/mihir/Documents/Summer/Models/combined_without_cpd43/ficks.xlsx"

active = macrophage
######################################
######################################

######################################
# EXCEL SHEET PARSED HERE  #
######################################
reactions = pd.read_excel(active, sheet_name = 1, skiprows = 1, header = 0)
species = pd.read_excel(active, sheet_name = 0, skiprows = 1, header = 0)
pmid = reactions["PMID"].tolist()
reactions = reactions[["Rule", "Weight", "n", "EC50"]]
species = species[["ID", "Yinit", "Ymax", "tau"]]
node_ID = species["ID"].tolist()
Yinit = species["Yinit"].tolist()
Ymax = species["Ymax"].tolist()
Ymax2 = species["Ymax"].tolist()
tau = species["tau"].tolist()

######################################
###### SIMULATOR FUNCTIONS HERE ######
######################################
# this function splits each reaction into it"s individual reactors, and returns an array with each reaction"s
# reactors in an array (I believe)
def get_reactors(reaction):
    reac_split = reaction.split(" ")
    reactors = []
    for k in reac_split:
        if k != "&" and k != "=>":
            reactors.append(k)
    return reactors[:-1]

# returns both reactors, unlike above function because Ficks requires the concentration of both species
def get_reactors_for_ficks(reaction):
    reac_split = reaction.split(" ")
    reactors = []
    for k in reac_split:
        if k != "&" and k != "=>" and k != "===>":
            reactors.append(k)
    return reactors

# activation function as defined by equations in PMID: 21087478
def Hill(reactor, n, EC50):
    # equation 1.3
    B = (EC50**n-1)/(2*EC50**n-1)
    C = (B-1)**(1/n)
        # if the first reactor has the prefix of ! then it is an inhibition reaction
        # equation 1.2
    if reactor[0] == "!":
        return (1-B*globals()["{}".format(reactor[1:])]**n/(C**n + globals()["{}".format(reactor[1:])]**n))
    else:
        return B*globals()["{}".format(reactor)]**n/(C**n + globals()["{}".format(reactor)]**n)


# returns the rate of change utilizing Ficks Second Law equation
def Ficks(transIn, transOut):
    # get C1 to be number of molecules of first species, C2 is number of molecule of second species in reaction
    conc_one = globals()["{}".format(transIn)]
    conc_two = globals()["{}".format(transOut)]
    deltaC = conc_two - conc_one
    dCoef = -0.05
    rate = dCoef * deltaC
    globals()["{}".format(transIn)] = globals()["{}".format(transIn)] - rate
    return rate

# if there are multiple possibilities for activation, this function is used
def OR(reaction_list):
    tera = (-1)**(len(reaction_list)+1)
    for k in reaction_list:
        weight, n, EC50 = reaction_list[k]
        final = weight
        for j in get_reactors(k):
            final *= Hill(j, n, EC50)
        tera *= (final-1)
    tera +=1
    return tera

# returns derivate (rate of change) for each species at each time point
def inte(state, t, reaction_dict):
    # setter for the state of each node
    for i in range(len(node_ID)):
        globals()["{}".format(node_ID[i])] = state[i]
    # for every node in the reaction
    for i in range(len(node_ID)):
        # TF represents a base reaction rate (I believe)
        TF = 1
        # if there is only one possible reaction
        if len(reaction_dict[node_ID[i]]) == 1:
            # get reaction string
            single_reaction = list(reaction_dict[node_ID[i]].keys())[0]
            # check reaction string and run if it is a Ficks diffusion reaction
            if "===>" in single_reaction:
                reactors = get_reactors_for_ficks(single_reaction)
                rate = Ficks(reactors[0], reactors[1])
                globals()["{}".format(node_ID[i] + "d")] = rate
            else:
                TF = 1
                # create a list of reactors from the reaction dictonary
                reactors = get_reactors(list(reaction_dict[node_ID[i]].keys())[0])
                weight, n, EC50 = reaction_dict[node_ID[i]][list(reaction_dict[node_ID[i]].keys())[0]]
                # for each inital reactor, get activation function from Hill
                for j in reactors:
                    TF *= Hill(j, n, EC50)
                # assignment to derivative of rate of change of node_ID[i]
                # equation derived from 1.1 from PMID: 21087478
                globals()["{}".format(node_ID[i] + "d")] = (TF*weight*Ymax[i]-globals()["{}".format(node_ID[i])])/tau[i]
        # otherwise, there are two possible reactions
        else:
            TF = OR(reaction_dict[node_ID[i]])
            globals()["{}".format(node_ID[i] + "d")] = (TF*Ymax[i]-globals()["{}".format(node_ID[i])])/tau[i]
    return [globals()["{}".format(k + "d")] for k in node_ID]

def hill_simulation(t, state0, reaction_dict):
    # state0 gets passed into inte() as the value of each species
    # odeint calls inte() for as many timepoints as t specifies
    yHill_ss = odeint(inte, state0, t, args = (reaction_dict,))
    print("Hill Finished\n")
    return yHill_ss

###################-------------------------------###################
###################----------Loading Data---------###################
###################-------------------------------###################

reaction_dict = collections.defaultdict(dict)
for k in range(len(reactions)):
    node = reactions.loc[k, "Rule"].split(" ")
    reaction_dict[node[-1]][reactions.loc[k, "Rule"]] = reactions.loc[k, ["Weight", "n", "EC50"]].tolist()


species_dict = dict()
for k in range(len(species)):
    #lis = species.loc[k, ["Yinit", "Ymax", "tau"]].tolist()
    species_dict[species.loc[k, "ID"]] = species.loc[k, ["Yinit", "Ymax", "tau"]].tolist()

state0 = []
for k in range(len(node_ID)):
    state0.append(Yinit[k])  #solve_ivp

######################################
# SET DISPLAY/EXPORT PARAMETERS HERE #
######################################
t = np.arange(0.0, 100, 0.01)
yHill_ss = hill_simulation(t, state0, reaction_dict)
whatToDisplay = 0
whatToDisplayTwo = 49
whatToExport = [10, 24, 25, 48, 49, 50, 55, 56, 57, 115, 116, 117]
timepointToExport = 1000
exportDataLocation = "../data/" + str(ntpath.basename(str(os.path.splitext(active)[0]))) + "_alldata.csv"
exportSensitivityAnalysisDataLocation = "../data/sensitivity_analysis/"
    + str(ntpath.basename(str(os.path.splitext(active)[0]))) + "_raw_sa_data.csv"
knockdownPercentage = 0.1
######################################

# number of timepoints to display
k = 10000

# Code to display graph
def displayGraph(indexOne, indexTwo, simData):
    plt.figure(figsize=(12,4))
    plt.subplot(121)
    plt.plot(t[:k], yHill_ss[:k,indexOne], label = node_ID[indexOne])
    plt.legend(loc="best")
    plt.subplot(122)
    plt.plot(t[:k], yHill_ss[:k,indexTwo], label = node_ID[indexTwo])
    plt.legend(loc="best")
    plt.show()

# Code to export a single species as a CSV
def exportSingleSpecies(whatToExport, simData):
    for eachSpecies in whatToExport:
        csvTitle = ("data/"+ node_ID[eachSpecies] + "_with_cpd43.csv")
        headerTitle = ("time," + node_ID[eachSpecies])
        data = np.transpose([t[:k], simData[:k, eachSpecies]])
        np.savetxt(csvTitle, data, delimiter=",", header=headerTitle)

# Code to export all data to a CSV
def exportAllData(exportLocation, simData):
    csv = open(exportLocation, "w")
    columnTitleRow = "time, "
    for species in node_ID:
        columnTitleRow += species + ","
    csv.write(columnTitleRow + "\n")
    timepoint_num = 0
    for timepoint in t.astype(str):
        csv.write(timepoint + ",")
        for species in range(len(node_ID)):
            csv.write(simData[timepoint_num,species].astype(str) + ",")
        timepoint_num += 1
        csv.write("\n")

# Code to export a given time point of a simulation
def exportLastTimePoint(exportLocation, simData, timepoint):
    csv = open(exportLocation, "w")
    csv.write("species, data")
    csv.write("\n")
    for species in range(len(node_ID)):
        csv.write(node_ID[species] + ",")
        csv.write(simData[timepoint,species].astype(str) + ",")
        csv.write("\n")

# The following functions are written as helper functions for the runAutoSensitivity function
# Code to write an initial control to the specified export location
def exportControlLastTimepoint(exportLocation, simData):
    csv = open(exportLocation, "w")
    csv.write("index, species, control")
    csv.write("\n")
    for species in range(len(node_ID)):
        csv.write(str(species) + "," + node_ID[species] + "," + simData[1000,species].astype(str) + ",")
        csv.write("\n")

# Code to write additional species' data to the csv
def exportForSensitivity(exportLocation, simData, knockdownSpecies):
    exportData = pd.read_csv(exportLocation)
    exportData.set_index('index', inplace=True)
    appendList = addValuesToList(simData)
    exportData[str(node_ID[knockdownSpecies])] = appendList
    exportData.to_csv(exportLocation)

# Code to append simulation data to a list so that it can be written to the pandas DataFrame
def addValuesToList(simData):
    list = []
    for species in range(len(node_ID)):
        list.append(simData[1000,species])
    return list

# Code that runs hill simulations with each Ymax knocked down to user-specified parameter, with an initial control
# simulation with no parameters changed
def runAutoSensitivity(knockdownPercentage, saLocation):
    controlSim = hill_simulation(t, state0, reaction_dict)
    exportControlLastTimepoint(saLocation, controlSim)
    for species in range(len(node_ID)):
        originalYMax = Ymax[species]
        newYmax = originalYMax * knockdownPercentage
        Ymax[species] = originalYMax * knockdownPercentage
        kdData = hill_simulation(t, state0, reaction_dict)
        exportForSensitivity(saLocation, kdData, species)
        Ymax[species] = originalYMax

######################################
## DISPLAY/EXPORT FUNCS CALLED HERE ##
######################################
runAutoSensitivity(knockdownPercentage, exportSensitivityAnalysisDataLocation)
# exportLastTimePoint(exportDataLocation, yHill_ss, timepointToExport)
# exportSingleSpecies(whatToExport, yHill_ss)
# exportAllData(exportDataLocation, yHill_ss)
# displayGraph(whatToDisplay, whatToDisplayTwo, yHill_ss)
######################################
