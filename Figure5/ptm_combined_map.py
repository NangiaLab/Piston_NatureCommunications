#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 13:26:38 2021

@author: kmpiston
"""
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import math


plt.rcParams['mathtext.default'] = 'tt'
#General Plotting Properties
#sns.set_style('ticks')
#sns.set_context('notebook',font_scale=2)


#Combined Heatmap of properties


sys = input("System name: ")
#paths to DF to load
path = "/home/kmpiston/DNA_analy/HEATMAP_ALL/"+sys+"/"
c_path = "/home/kmpiston/DNA_analy/CONTACTS/resDFs/"+sys+"/"
h_path = '/home/kmpiston/DNA_analy/HBOND/resDFs/'+sys+"/"


#Custom Colorbar Creation
cmap = mpl.colors.ListedColormap(['white','royalblue','xkcd:saffron','xkcd:fresh green','xkcd:scarlet'])
bounds = [0,1,2,3,4,5]
norm= mpl.colors.BoundaryNorm(bounds, cmap.N)

#Load Dfs

dna_dfs = []
hbond_dfs = []
CL_dfs = []
NA_dfs = []

N = [*range(1,39)] # residue numbers

for n in N:
    dna_dfs.append(pd.read_pickle(c_path+'contacts'+str(n)+'.pkl'))
    hbond_dfs.append(pd.read_pickle(h_path+'hbond'+str(n)+'.pkl'))
    CL_dfs.append(pd.read_pickle(c_path+'ions_CL'+str(n)+'.pkl'))
    NA_dfs.append(pd.read_pickle(c_path+'ions_NA'+str(n)+'.pkl'))

dt = math.ceil(501/len(dna_dfs[0].columns))

no_c=[] # empty variable to store residues with no contacts
#change names column names 
for n in N:
    dnas = len(dna_dfs[n-1].index)
    inds = []
    if dnas == 1:
        v = dna_dfs[n-1].index[0]
        if v == 'Contacts':
            inds=[' ']
            no_c.append(n) #flag residues that have no contact so they arn't deleted later with the zero rows
        else:
            inds.append('d'+str(v))
    else:
        for d in range(dnas):
            v = dna_dfs[n-1].index[d]
            inds.append('d'+str(v))
    dna_dfs[n-1].index = inds
    hbond_dfs[n-1].index = ['w-hb']
    CL_dfs[n-1].index = ['Cl$^-$']
    NA_dfs[n-1].index = ['Na$^+$']

#Path until new trajectory is rendered for ACT, right now the traj has 267 frames and it should have only 251, we will delete the first 16 columns
cdrop = [*range(0,16)]

if sys == 'act':
    for n in N:
        dna_dfs[n-1].drop(dna_dfs[n-1].columns[cdrop], axis=1, inplace=True)
        hbond_dfs[n-1].drop(hbond_dfs[n-1].columns[cdrop], axis=1, inplace=True)
        CL_dfs[n-1].drop(CL_dfs[n-1].columns[cdrop], axis=1, inplace=True)
        NA_dfs[n-1].drop(NA_dfs[n-1].columns[cdrop], axis=1, inplace=True)
    
#adjust frequency of data if dt != 5 
if dt != 10:
    for n in N:
        dna_dfs[n-1] = dna_dfs[n-1][dna_dfs[n-1].columns[::5]]
        hbond_dfs[n-1] = hbond_dfs[n-1][hbond_dfs[n-1].columns[::5]]
        CL_dfs[n-1] = CL_dfs[n-1][CL_dfs[n-1].columns[::5]]
        NA_dfs[n-1] = NA_dfs[n-1][NA_dfs[n-1].columns[::5]]
    dt = 10

times = [*range(0,501,dt)]

#combine dfs to form one heatmap and drop any DNA df columns with all Zeros
for n in N:
    dna_dfs[n-1].columns=times 
    hbond_dfs[n-1].columns=times
    CL_dfs[n-1].columns = times
    NA_dfs[n-1].columns = times
#Drop columns with all zeros except those that have no contacts 

#get list of values that have contacts by comparing N  and no_c

z = list(set(N).symmetric_difference(set(no_c)))
for n in z:
    if dna_dfs[n-1].max().max() > 0:
        dna_dfs[n-1] = dna_dfs[n-1][(dna_dfs[n-1].T !=0).any()] # delete zeros 
    else: 
        if len(dna_dfs[n-1].index) >1:
            num = len(dna_dfs[n-1].index) - 1
            while num > 0:
                dna_dfs[n-1].drop(index=dna_dfs[n-1].index[-1], axis=0, inplace=True) #drop last row, this is in case of multiple empty rows
                num -= 1
        dna_dfs[n-1].index = [' ']

#Combine into one heatmap and adjust values for colors and store size of each residue
sizes = []

for n in N:
    df = pd.concat([dna_dfs[n-1],hbond_dfs[n-1],CL_dfs[n-1],NA_dfs[n-1]])
    sizes.append(len(df.index))
        
   
    
    #change df to binary to adjust values for heatmap
    df = df.ge(1.0).astype(int) #threshold for up or down was set for 1.0
    
    
    
    #number of dna contact rows
    dna_rows = (len(df.index)-3)
    dna_rows
    
    d = df.iloc[0:dna_rows,:].astype(float)
    for v in d.values:
        v *= .2
    df.iloc[0:dna_rows,:] = d
    
    #Hbond color
    h = df.iloc[[dna_rows],:].astype(float)
    for v in h.values:
        v *=.4
    df.iloc[[dna_rows],:]= h
    
    #CL Color Green
    c = df.iloc[[dna_rows+1],:].astype(float)
    for v in c.values:
        v *=.6
    df.iloc[[dna_rows+1],:]= c
    
    #NA Color Green
    Na = df.iloc[[dna_rows+2],:].astype(float)
    for v in Na.values:
        v *=.8
    df.iloc[[dna_rows+2],:]= Na
    
    
    fig_dims= (10,3)
    plt.rcParams['font.size'] = 12
    fig, ax = plt.subplots(figsize = fig_dims, edgecolor = 'k')
    sns.heatmap(df, 
                cmap=cmap,
                vmin=0,
                vmax=.8,
                cbar=False,
                linewidths=1.5, 
                linecolor='white', 
                square=False)
    plt.xticks(rotation = 90,fontsize=18, fontname='Arial', fontstretch='condensed', weight='bold')
    plt.yticks(rotation = 0, fontsize=18, fontname= 'Arial',fontstretch='condensed', weight = 'bold')
    #ax.set_xticklabels(labels= ax.get_xticklabels(), ha = 'left')
    ax.tick_params(direction='in', left=True, top = True, right = True, width=2, length =10)
    h = ax.get_yticks()
    l = len(h)
    w = ax.get_xticks()
    ax.hlines((h[l-4]+h[l-3])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray") #gridlines
    ax.hlines((h[l-3]+h[l-2])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray")
    ax.hlines((h[l-2]+h[l-1])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray")
    ax.spines["top"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["top"].set_linewidth(5)
    ax.spines["bottom"].set_linewidth(5)
    ax.spines["left"].set_linewidth(5)
    ax.spines["right"].set_linewidth(5)
    #plt.xlim((0,500))
    plt.tight_layout() 
    plt.savefig(path+str(sys)+str(n)+".png", dpi=300, transparent=True)
    df.to_pickle(path+str(n)+'.pkl')

#Combine to one big heatmap
df_whole = pd.read_pickle(path+str(1)+'.pkl')

N2 = [*range(2,39)] # residue numbers 2 to end

for r in N2:
    df2 = pd.read_pickle(path+str(r)+'.pkl')
    df_whole = pd.concat([df_whole,df2])

#remove last value of sizes for loop so border can be drawn
sizes.pop()

#Plot Whole Fig
fig_dims= (30,50)
plt.rcParams['font.size'] = 18
fig, ax = plt.subplots(figsize = fig_dims, edgecolor = 'k')
sns.heatmap(df_whole, 
            cmap=cmap,
            vmin=0,
            vmax=.8,
            cbar=False,
            linewidths=1.5, 
            linecolor='white', 
            square=False)
y = df_whole.index.to_list() # setting y axis to every value
ax.set_yticks((np.arange(len(y))+.5))
ax.set_yticklabels(y)
plt.xticks(rotation = 90,fontsize=18, fontname='Arial', fontstretch='condensed', weight='bold')
plt.yticks(rotation = 0, fontsize=18, fontname= 'Arial',fontstretch='condensed', weight = 'bold')
#ax.set_xticklabels(labels= ax.get_xticklabels(), ha = 'left')
ax.tick_params(direction='in', left=True, top = True, right = True, width=2, length =10)
h = ax.get_yticks()
l = len(h)
w = ax.get_xticks()
flag = 0 #Flag what row we are moving through to add the gridline
for s in sizes:
    flag = flag + s
    ax.hlines((h[flag]+h[flag-1])/2,w[-1]+.5,w[0]-.5, linestyle='solid', linewidth=2, color="black") #gridlines
    ax.hlines((h[flag-4]+h[flag-3])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray") #gridlines
    ax.hlines((h[flag-3]+h[flag-2])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray")
    ax.hlines((h[flag-2]+h[flag-1])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray")
#Add last row of lines
ax.hlines((h[l-4]+h[l-3])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray") #gridlines
ax.hlines((h[l-3]+h[l-2])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray")
ax.hlines((h[l-2]+h[l-1])/2,w[-1]+.5,w[0]-.5, linestyle='dotted', linewidth=2, color="dimgray")
#ax.axes.get_yaxis().set_ticks([])

ax.spines["top"].set_visible(True)
ax.spines["bottom"].set_visible(True)
ax.spines["left"].set_visible(True)
ax.spines["right"].set_visible(True)
ax.spines["top"].set_linewidth(5)
ax.spines["bottom"].set_linewidth(5)
ax.spines["left"].set_linewidth(5)
ax.spines["right"].set_linewidth(5)
plt.tight_layout() 
plt.savefig(path+str(sys)+"_whole.png", dpi=300, transparent=True)
df_whole.to_pickle(path+ "whole.pkl")
