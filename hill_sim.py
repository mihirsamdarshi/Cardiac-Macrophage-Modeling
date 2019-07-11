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
from scipy.optimize import minimize

############################
# SET THE EXCEL SHEET HERE #
############################
macrophage = '/Users/mihir/Documents/Summer/Models/macrophage_model.xlsx'
no_inhibition = '/Users/mihir/Documents/Summer/Models/macrophage_model_no_inhibition.xlsx'
step_by_step = '/Users/mihir/Documents/Summer/Models/step_by_step_model.xlsx'
og_fibroblast = '/Users/mihir/Documents/Summer/Models/original_fibroblast_model.xls'
og_cardiomyocyte = '/Users/mihir/Documents/Summer/Models/original_cardiomyocyte_model.xls'

active = step_by_step
############################
############################

############################
# EXCEL SHEET PARSED HERE  #
############################
reactions = pd.read_excel(active, sheet_name = 1, skiprows = 1, header = 0)
species = pd.read_excel(active, sheet_name = 0, skiprows = 1, header = 0)
pmid = reactions['PMID'].tolist()
reactions = reactions[['Rule', 'Weight', 'n', 'EC50']]
species = species[['ID', 'Yinit', 'Ymax', 'tau']]
node_ID = species['ID'].tolist()
Yinit = species['Yinit'].tolist()
Ymax = species['Ymax'].tolist()
Ymax2 = species['Ymax'].tolist()
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

# read initial state
state0 = []
for k in range(len(node_ID)):
    state0.append(Yinit[k])  #solve_ivp

# Set Time points here
t = np.arange(0.0, 10, 0.1)

############################
# SIMULATOR FUNCTIONS HERE #
############################
# this function splits each reaction into it's individual reactors, and returns an array with each reaction's
# reactors in an array (I believe)
def get_reactors(reaction):
    reac_split = reaction.split(' ')
    reactors = []
    for k in reac_split:
        if k != '&' and k != '=>' and k != '===>':
            reactors.append(k)
    return reactors[:-1]

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
    K = (B-1)**(1/n)
    # if the first reactor has the prefix of ! then it is an inhibition reaction
    # equation 1.2
    if reactor[0] == '!':
        return (1-B*globals()['{}'.format(reactor[1:])]**n/(K**n + globals()['{}'.format(reactor[1:])]**n))
    else:
        return B*globals()['{}'.format(reactor)]**n/(K**n + globals()['{}'.format(reactor)]**n)

def Ficks(transIn, transOut):
    # get C1 to be number of molecules of first species, C2 is number of molecule of second species in reaction
    conc_one = globals()['{}'.format(transIn)]
    conc_two = globals()['{}'.format(transOut)]
    deltaC = conc_two - conc_one
    deltaX = 1
    D = 1

    return -D * deltaC / deltaX

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
            if any('===>' in string for string in single_reaction):
                reactors = get_reactors_for_ficks(single_reaction)
                rate = Ficks(reactors[0], reactors[1])
                globals()['{}'.format(node_ID[i] + 'd')] = rate
                print('true')
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
    print('Hill Finished')
    return yHill_ss

yHill_ss = hill_simulation(t, state0, reaction_dict)

############################
# SET THE EXCEL SHEET HERE #
############################
whatToDisplay = 3
############################
############################

# Code to display graph
k = 12000
plt.figure(figsize=(12,4))
plt.subplot(121)
plt.plot(t[:k], yHill_ss[:k,whatToDisplay], label = node_ID[whatToDisplay])
plt.legend(loc='best')
plt.show()
