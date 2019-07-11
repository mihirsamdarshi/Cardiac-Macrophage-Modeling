#!/usr/bin/env python
# coding: utf-8

# ## Cardiac Regulatory Network for CMRG
# ### CRNC developed by Shulin Cao
# ### Only for use in CMRG

# In[26]:


###------Network Simulator------###
###------Shulin Cao------###
###------CMRG, UC San Diego------###

###import packages###
import pandas as pd
import collections
import timeit
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.mlab as mlab
import statistics
import random
import numpy.linalg
import numpy as np
import sys
from scipy.optimize import minimize
elapsed_time = timeit.default_timer()
from sklearn.linear_model import LinearRegression
from sklearn import cluster
import seaborn as sns
sns.set()
from sklearn import datasets
from sklearn.metrics import r2_score
import csv

############################
# SET THE EXCEL SHEET HERE #
############################
regular = '/Users/mihir/Documents/Summer/Models/macrophage_model.xlsx'
no_inhibition = '/Users/mihir/Documents/Summer/Models/macrophage_model_no_inhibition.xlsx'
stepbystep = '/Users/mihir/Documents/Summer/Models/stepbystep.xlsx'

active = stepbystep
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

reaction_dict = collections.defaultdict(dict)
for k in range(len(reactions)):
    node = reactions.loc[k, 'Rule'].split(' ')
    reaction_dict[node[-1]][reactions.loc[k, 'Rule']] = reactions.loc[k, ['Weight', 'n', 'EC50']].tolist()

species_dict = dict()
for k in range(len(species)):
    species_dict[species.loc[k, 'ID']] = species.loc[k, ['Yinit', 'Ymax', 'tau']].tolist()

geneinput = {i:0 for i in node_ID[98:]}
for j in node_ID[98:]:
    for k in reaction_dict[j]:
        if '!' in k:
            geneinput[j] = 0.5
            break

# read initial state
state0 = []
for k in range(len(node_ID)):
    state0.append(Yinit[k])  #solve_ivp

# Set Time points here
t = np.arange(0, 4, 1)

############################
# SIMULATOR FUNCTIONS HERE #
############################
# this function splits each reaction into it's individual reactos, and returns an array with each reaction's
# reactors in an array (I believe)
def get_reactors(reac):
    reac_split = reac.split(' ')
    reactors = []
    for k in reac_split:
        if k != '&' and k != '=>' and k != '===>':
            reactors.append(k)
    return reactors[:-1]

def Hill(reactor, n, EC50):
    B = (EC50**n-1)/(2*EC50**n-1)
    C = (B-1)**(1/n)
    # if the first reactor has the prefix of ! then it is an inhibition reaction
    if reactor[0] == '!':
        return (1-B*globals()['{}'.format(reactor[1:])]**n/(C**n + globals()['{}'.format(reactor[1:])]**n))
    else:
        return B*globals()['{}'.format(reactor)]**n/(C**n + globals()['{}'.format(reactor)]**n)

# def Ficks(reaction_dict[node_ID[i]]):
#     Assuming c and x are numpy arrays of equal size and D is a scalar
#     get C1 to be number of molecules of first species, C2 is number of molecule of second species in reaction
#     deltaC = np.diff(c)
#     deltaX = distance between macrophage and fibroblast
#
#   return -D * deltaC / deltaX

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
        # create a list of reactions
        allReactions = list(reaction_dict[node_ID[i]].keys())[0]
        # if this is a Ficks diffusion reaction
        if any('===>' in string for string in reactors):
            print('true')
        # else if there is only one possible reaction
        elif len(reaction_dict[node_ID[i]]) == 1:
            # create a list of reactors from the reaction dictonary
            reactors = get_reactors(list(allReactions)
            weight, n, EC50 = reaction_dict[node_ID[i]][allReactions]
            TF = 1
            for j in reactors:
                TF *= Hill(j, n, EC50)
            globals()['{}'.format(node_ID[i] + 'd')] = (TF*weight*Ymax[i]-globals()['{}'.format(node_ID[i])])/tau[i]
        else:
            TF = OR(reaction_dict[node_ID[i]])
            globals()['{}'.format(node_ID[i] + 'd')] = (TF*Ymax[i]-globals()['{}'.format(node_ID[i])])/tau[i]
    # print(list(k for k in node_ID))
    # print(list([globals()['{}'.format(k + 'd')] for k in node_ID]))
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
whatToDisplay = 4
############################
############################

# k = 12000
# plt.figure(figsize=(12,4))
# plt.subplot(121)
# plt.plot(t[:k], yHill_ss[:k,whatToDisplay], label = node_ID[whatToDisplay])
# plt.legend(loc='best')
# plt.show()
