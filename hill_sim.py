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
sns.set()

######################################
# SET THE EXCEL SHEET HERE #
######################################
macrophage = '/Users/mihir/Documents/Summer/Models/macrophage_model.xlsx'
no_inhibition = '/Users/mihir/Documents/Summer/Models/macrophage_model_no_inhibition.xlsx'
step_by_step = '/Users/mihir/Documents/Summer/Models/step_by_step_model.xlsx'
og_fibroblast = '/Users/mihir/Documents/Summer/Models/original_fibroblast_model.xls'
og_cardiomyocyte = '/Users/mihir/Documents/Summer/Models/original_cardiomyocyte_model.xlsx'
ficks = '/Users/mihir/Documents/Summer/Models/ficks.xlsx'
combined = '/Users/mihir/Documents/Summer/Models/combined_model.xlsx'

active = ficks
######################################
######################################

######################################
# EXCEL SHEET PARSED HERE  #
######################################
reactions = pd.read_excel(active, sheet_name = 1, skiprows = 1, header = 0)
species = pd.read_excel(active, sheet_name = 0, skiprows = 1, header = 0)
pmid = reactions['PMID'].tolist()
reactions = reactions[['Rule', 'Weight', 'n', 'EC50']]
species = species[['ID', 'Yinit', 'Ymax', 'tau']]
node_ID = species['ID'].tolist()
Yinit = species['Yinit'].tolist()
Ymax = species['Ymax'].tolist()
tau = species['tau'].tolist()

# create a dictionary of all the reactions
reaction_dict = collections.defaultdict(dict)
for k in range(len(reactions)):
    node = reactions.loc[k, 'Rule'].split(' ')
    reaction_dict[node[-1]][reactions.loc[k, 'Rule']] = reactions.loc[k, ['Weight', 'n', 'EC50']].tolist()

# create a dictionary of all the species
species_dict = dict()
for k in range(len(species)):
    species_dict[species.loc[k, 'ID']] = species.loc[k, ['Yinit', 'Ymax', 'tau']].tolist()

# read and set the initial state based on Yinit from Excel sheet
state0 = []
for k in range(len(node_ID)):
    state0.append(Yinit[k])  #solve_ivp
############################
############################

######################################
###### SIMULATOR FUNCTIONS HERE ######
######################################
# this function splits each reaction into it's individual reactors, and returns an array with each reaction's
# reactors in an array (I believe)
def get_reactors(reaction):
    reac_split = reaction.split(' ')
    reactors = []
    for k in reac_split:
        if k != '&' and k != '=>':
            reactors.append(k)
    return reactors[:-1]

# returns both reactors, unlike above function because Ficks requires the concentration of both species
def get_reactors_for_ficks(reaction):
    reac_split = reaction.split(' ')
    reactors = []
    for k in reac_split:
        if k != '&' and k != '=>' and k != '===>':
            reactors.append(k)
    return reactors

# activation function as defined by equations in PMID: 21087478
def Hill(reactor, n, EC50, reaction_flux_values):
    # equation 1.3
    B = (EC50**n-1)/(2*EC50**n-1)
    C = (B-1)**(1/n)
        # if the first reactor has the prefix of ! then it is an inhibition reaction
        # equation 1.2
    if reactor[0] == '!':
        reactor_value = reaction_flux_values[reactor[1:]]
        return (1-B*reactor_value**n)/ (C**n + reactor_value**n)
    else:
        reactor_value = reaction_flux_values[reactor]
        return (B*reactor_value**n)/ (C**n + reactor_value**n)


# returns the rate of change utilizing Ficks Second Law equation
def Ficks(transIn, transOut, reaction_flux_values):
    # get C1 to be number of molecules of first species, C2 is number of molecule of second species in reaction
    conc_one = reaction_flux_values[transIn]
    print(transIn + ": " + str(conc_one))
    conc_two = reaction_flux_values[transOut]
    print(transOut + ": " + str(conc_two))
    deltaC = conc_two - conc_one
    print('deltaC: ' + str(deltaC))
    dCoef = -0.05
    rate = dCoef * deltaC
    print('rate: ' + str(rate))
    reaction_flux_values[transIn] = reaction_flux_values[transIn] - rate
    return rate

# if there are multiple possibilities for activation, this function is used
def OR(reaction_list, reaction_flux_values):
    tera = (-1)**(len(reaction_list)+1)
    for k in reaction_list:
        weight, n, EC50 = reaction_list[k]
        final = weight
        for j in get_reactors(k):
            final *= Hill(j, n, EC50, reaction_flux_values)
        tera *= (final-1)
    tera +=1
    return tera

# returns derivate (rate of change) for each species at each time point
def inte(state, t, reaction_dict):
    # setter for the state of each node
    reaction_flux_values = dict()
    for i in range(len(node_ID)):
        reaction_flux_values.update([(node_ID[i], state[i])])
    print(reaction_flux_values)
    # for every node in the reaction
    for i in range(len(node_ID)):
        # TF represents a base reaction rate (I believe)
        TF = 1
        print(reaction_flux_values)
        # if there is only one possible reaction
        if len(reaction_dict[node_ID[i]]) == 1:
            # get reaction string
            single_reaction = list(reaction_dict[node_ID[i]].keys())[0]
            # check reaction string and run if it is a Ficks diffusion reaction
            if '===>' in single_reaction:
                print('Ficks')
                reactors = get_reactors_for_ficks(single_reaction)
                rate = Ficks(reactors[0], reactors[1], reaction_flux_values)
                print('rate of ' + single_reaction + ': ' + str(rate))
                current_value = reaction_flux_values[node_ID[i]]
                print('current value of ' + node_ID[i] + ' : ' + str(current_value))
                current_value += rate
                print('after addition ' + node_ID[i] + ' : ' + str(current_value))
                reaction_flux_values[node_ID[i]] = current_value
            else:
                print('Hill')
                # create a list of reactors from the reaction dictonary
                reactors = get_reactors(list(reaction_dict[node_ID[i]].keys())[0])
                print(reactors)
                weight, n, EC50 = reaction_dict[node_ID[i]][list(reaction_dict[node_ID[i]].keys())[0]]
                # for each inital reactor, get activation function from Hill
                for j in reactors:
                    TF *= Hill(j, n, EC50, reaction_flux_values)
                # assignment to derivative of rate of change of node_ID[i]
                # equation derived from 1.1 from PMID: 21087478
                current_value = reaction_flux_values[node_ID[i]]
                reaction_flux_values[node_ID[i]] = (TF*weight*Ymax[i]-current_value)/tau[i]
        # otherwise, there are two possible reactions
        else:
            print('hello')
            current_value = reaction_flux_values[node_ID[i]]
            TF = OR(reaction_dict[node_ID[i]], reaction_flux_values)
            reaction_flux_values[node_ID[i]] = (TF*weight*Ymax[i]-current_value)/tau[i]
    print('\n')
    return [reaction_flux_values[k] for k in node_ID]

def hill_simulation(t, state0, reaction_dict):
    # state0 gets passed into inte() as the value of each species
    # odeint calls inte() for as many timepoints as t specifies
    yHill_ss = odeint(inte, state0, t, args = (reaction_dict,))
    print('Hill Finished\n')
    return yHill_ss

######################################
# SET DISPLAY/EXPORT PARAMETERS HERE #
######################################
t = np.arange(0.0, 10, 1)
yHill_ss = hill_simulation(t, state0, reaction_dict)
whatToDisplay = 0
whatToDisplayTwo = 1
whatToExport = 0
exportDataLocation = "data/allData.csv"
knockdownPercentage = 0.5
######################################

# number of timepoints to display
k = 6000

# Code to display graph
def displayGraph(whatToDisplay, simData):
    plt.figure(figsize=(12,4))
    plt.subplot(121)
    plt.plot(t[:k], simData[:k,whatToDisplay], label = node_ID[whatToDisplay])
    plt.legend(loc='best')
    plt.show()

# Code to export a single species as a CSV
def exportSingleSpecies(whatToExport, simData):
    csvTitle = ("Data/"+ node_ID[whatToExport] + "_inhibited.csv")
    headerTitle = ('time,' + node_ID[whatToExport])
    data = np.transpose([t[:k], simData[:k,whatToExport]])
    np.savetxt(csvTitle, data, delimiter=",", header=headerTitle)

# Code to export all data to a CSV
def exportAllData(exportLocation, simData):
    csv = open(exportLocation, "w")
    columnTitleRow = "time, "
    for species in node_ID:
        columnTitleRow += species + ","
    csv.write(columnTitleRow + '\n')
    timepoint_num = 0
    for timepoint in t.astype(str):
        csv.write(timepoint + ',')
        for species in range(len(node_ID)):
            csv.write(simData[timepoint_num,species].astype(str) + ",")
        timepoint_num += 1
        csv.write('\n')

# Code that runs hill simulations with each Ymax knocked down to user-specified parameter
def runAutoSensitivity(knockdownPercentage):
    for species in range(len(node_ID)):
        originalYMax = Ymax[species]
        newYmax = originalYMax * knockdownPercentage
        Ymax[species] = originalYMax * knockdownPercentage
        kdData = hill_simulation(t, state0, reaction_dict)
        saLocation = "data/sensitivity_analysis/sa_" + str(knockdownPercentage) + "_" +  node_ID[species] + ".csv"
        exportAllData(saLocation, kdData)
        Ymax[species] = originalYMax

######################################
## DISPLAY/EXPORT FUNCS CALLED HERE ##
######################################
# runAutoSensitivity(knockdownPercentage)
# exportSingleSpecies(whatToExport, yHill_ss)
# exportAllData(exportDataLocation, yHill_ss)
# displayGraph(whatToDisplay, yHill_ss)
# displayGraph(whatToDisplayTwo, yHill_ss)
######################################
######################################
