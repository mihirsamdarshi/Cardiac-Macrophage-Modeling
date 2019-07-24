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
np.seterr(all='warn')

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

active = combined
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
def Hill(reactor, n, EC50):
    # equation 1.3
    B = (EC50**n-1)/(2*EC50**n-1)
    C = (B-1)**(1/n)
        # if the first reactor has the prefix of ! then it is an inhibition reaction
        # equation 1.2
    if reactor[0] == '!':
        return (1-B*globals()['{}'.format(reactor[1:])]**n/(C**n + globals()['{}'.format(reactor[1:])]**n))
    else:
        return B*globals()['{}'.format(reactor)]**n/(C**n + globals()['{}'.format(reactor)]**n)


# returns the rate of change utilizing Ficks Second Law equation
def Ficks(transIn, transOut):
    # get C1 to be number of molecules of first species, C2 is number of molecule of second species in reaction
    conc_one = globals()['{}'.format(transIn)]
    conc_two = globals()['{}'.format(transOut)]
    deltaC = conc_two - conc_one
    deltaX = 1
    area = 1
    dCoef = 0.5
    rate = dCoef * -area * deltaC / deltaX
    globals()['{}'.format(transIn)] = globals()['{}'.format(transIn)] - rate
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
        globals()['{}'.format(node_ID[i])] = state[i]
    # for every node in the reaction
    for i in range(len(node_ID)):
        # TF represents a base reaction rate (I believe)
        TF = 1
        # if there is only one possible reaction
        if len(reaction_dict[node_ID[i]]) == 1:
            # get reaction string
            single_reaction = list(reaction_dict[node_ID[i]].keys())[0]
            # check reaction string and run if it is a Ficks diffusion reaction
            if '===>' in single_reaction:
                reactors = get_reactors_for_ficks(single_reaction)
                rate = Ficks(reactors[0], reactors[1])
                globals()['{}'.format(node_ID[i] + 'd')] = rate
            else:
                # create a list of reactors from the reaction dictonary
                reactors = get_reactors(list(reaction_dict[node_ID[i]].keys())[0])
                weight, n, EC50 = reaction_dict[node_ID[i]][list(reaction_dict[node_ID[i]].keys())[0]]
                # for each inital reactor, get activation function from Hill
                for j in reactors:
                    TF *= Hill(j, n, EC50)
                # assignment to derivative of rate of change of node_ID[i]
                # equation derived from 1.1 from PMID: 21087478
                globals()['{}'.format(node_ID[i] + 'd')] = (TF*weight*Ymax[i]-globals()['{}'.format(node_ID[i])])/tau[i]
        # otherwise, there are two possible reactions
        else:
            TF = OR(reaction_dict[node_ID[i]])
            globals()['{}'.format(node_ID[i] + 'd')] = (TF*Ymax[i]-globals()['{}'.format(node_ID[i])])/tau[i]
    return [globals()['{}'.format(k + 'd')] for k in node_ID]

def hill_simulation(t, state0, reaction_dict):
    # state0 gets passed into inte() as the value of each species
    # odeint calls inte() for as many timepoints as t specifies
    yHill_ss = odeint(inte, state0, t, args = (reaction_dict,))
    print('Hill Finished\n')
    return yHill_ss

######################################
# SET DISPLAY/EXPORT PARAMETERS HERE #
######################################
t = np.arange(0.0, 60, 0.1)
yHill_ss = hill_simulation(t, state0, reaction_dict)
whatToDisplay = 11
whatToDisplayTwo = 80
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
displayGraph(whatToDisplay, yHill_ss)
displayGraph(whatToDisplayTwo, yHill_ss)
######################################
######################################
