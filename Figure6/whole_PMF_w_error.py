#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 11:41:40 2021

@author: kmpiston
"""

#Paths to files

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy import interpolate
from scipy.signal import medfilt
import matplotlib.font_manager #fixes Arial font issue

sns.set_style('ticks')
sns.set_context('notebook',font_scale=3)

 

path =   ('/home/kmpiston/h3control_pull/bsResult.xvg') # Defining path
path_phos =  ('/home/kmpiston/Phosphorylation/charmm-gui-0130791649/gromacs/final_bsResult.xvg')
path_aly=  ('/home/kmpiston/H3K4_ALY/gromacs/final_bsResult.xvg')
path_mml =  ('/home/kmpiston/bridges_MML/gromacs/final_bsResult.xvg')
path_dml =  ('/home/kmpiston/DML_pull/charmm-gui-0993440136/gromacs/final_bsResult.xvg')
path_tml =  ('/home/kmpiston/TML/charmm-gui-0745037539/gromacs/final_bsResult.xvg')
path_act = ('/home/kmpiston/multi_pull/active/final_bsResult.xvg')
path_rep = ('/home/kmpiston/multi_pull/repressive/charmm-gui-2627332053/gromacs/final_bsResult.xvg')


#path_tml =('/home/kmpiston/TML/charmm-gui-0745037539/gromacs/50profile.xvg')

#Load
nm, e, yerr    = np.loadtxt(path, unpack=True)
nm2, e2, yerr2 = np.loadtxt(path_phos, unpack=True)
nm3, e3, yerr3 = np.loadtxt(path_aly, unpack=True)
nm4, e4, yerr4 = np.loadtxt(path_mml, unpack=True)
nm5, e5, yerr5 = np.loadtxt(path_dml, unpack=True)
nm6, e6, yerr6 = np.loadtxt(path_tml, unpack=True)
nm7, e7, yerr7 = np.loadtxt(path_act, unpack=True)
nm8, e8, yerr8 = np.loadtxt(path_rep, unpack=True)

#shift 
E1 = e - e.max()
E2 = e2 - e2.max()
E3 = e3 - e3.max()
E4 = e4 - e4.max()
E5 = e5 - e5.max()
E6 = e6 - e6.max()
E7 = e7 - e7.max()
E8 = e8 - e8.max()


#Shift Align and Smooth
x = [nm, nm2,nm3, nm4, nm5, nm6, nm7, nm8]
y = [E1, E2, E3, E4, E5, E6, E7, E8]

#Shift all
x_shift = []
for p in x:
    shift = p.min() - nm.min()
    x_shift.append(p-shift)

#Median filter
y_smooth = []
filter_length = 3 #FILTER TO MANIPULATE
for e in y:
    y_smooth.append(medfilt(e, filter_length))



nm = x_shift
E  = y_smooth
#Plot
fig,ax = plt.subplots(figsize=(20,15))
#plt.rcParams['axes.linewidth'] = 1 #set the value globally

#control
line_c,= ax.plot(nm[0],E[0], color ="black", linewidth = 6)
ax.fill_between(nm[0], E[0]-yerr, E[0]+yerr, color=line_c.get_color(), alpha =.3)
#Phosphorylation
line_c,= ax.plot(nm[1],E[1], color ="blue", linewidth = 6)
ax.fill_between(nm[1], E[1]-yerr2, E[1]+yerr2, color=line_c.get_color(), alpha =.3)
##Acetylation K4
line_c,= ax.plot(nm[2],E[2], color ="hotpink", linewidth = 6)
ax.fill_between(nm[2], E[2]-yerr3, E[2]+yerr3, color=line_c.get_color(), alpha =.3)
##Me1
line_c,= ax.plot(nm[3],E[3], color ="purple", linewidth = 6)
ax.fill_between(nm[3], E[3]-yerr4, E[3]+yerr4, color=line_c.get_color(), alpha =.3)
##Me2
line_c,= ax.plot(nm[4],E[4], color ="cyan", linewidth = 6)
ax.fill_between(nm[4], E[4]-yerr5, E[4]+yerr5, color=line_c.get_color(), alpha =.3)
##Me3
line_c,= ax.plot(nm[5],E[5], color ="orange", linewidth = 6)
ax.fill_between(nm[5], E[5]-yerr[5], E[5]+yerr6, color=line_c.get_color(), alpha =.3)
#ax.fill_between(nm6, E6-.35, E6+.35, color=line_c.get_color(), alpha =.3)

##4 Mes
line_c,= ax.plot(nm[6], E[6], color ="green", linewidth = 6)
ax.fill_between(nm[6], E[6]-yerr7, E[6]+yerr7, color=line_c.get_color(), alpha =.3)

line_c,= ax.plot(nm[7], E[7], color ="red", linewidth = 6)
ax.fill_between(nm[7], E[7]-yerr8, E[7]+yerr8, color=line_c.get_color(), alpha =.3)

#ax.legend(labels=["Control", "S10,28 Phosphorylation","K4 Acetylation", "K4 Me1", "K4 Me2", "K4 Me3", "4 Me"])
ax.tick_params(direction="in", width=6, length =24, labelsize = 48, pad=8, top=True, right=True)
plt.yticks(fontname='Arial')
ax.spines["top"].set_linewidth(10)
ax.spines["bottom"].set_linewidth(10)
ax.spines["left"].set_linewidth(10)
ax.spines["right"].set_linewidth(10)

#ax.set_tick_params(width=5)
#Despine
#sns.despine()


#plt.title("H3 Tail PMF Curves")
plt.xlabel("Discrete Coordinates (nm)", fontname='Arial')

plt.ylabel(' Potential ${\Delta}$E (kJ $mol^{-1}$)', fontname ='Arial')
#plt.savefig('/home/kmpiston/DNA_analy/PMF/new_pmf_combined.png', dpi=300, transparent=True)

#Create Insert Plot

fig,ax = plt.subplots(figsize=(20,20))
#plt.rcParams['axes.linewidth'] = 1 #set the value globally

#control
line_c,= ax.plot(nm[0],E[0], color ="black", linewidth = 10)
ax.fill_between(nm[0], E[0]-yerr, E[0]+yerr, color=line_c.get_color(), alpha =.3)
#Phosphorylation
line_c,= ax.plot(nm[1],E[1], color ="blue", linewidth = 10)
ax.fill_between(nm[1], E[1]-yerr2, E[1]+yerr2, color=line_c.get_color(), alpha =.3)
##Acetylation K4
line_c,= ax.plot(nm[2],E[2], color ="hotpink", linewidth = 10)
ax.fill_between(nm[2], E[2]-yerr3, E[2]+yerr3, color=line_c.get_color(), alpha =.3)
##Me1
line_c,= ax.plot(nm[3],E[3], color ="purple", linewidth = 10)
ax.fill_between(nm[3], E[3]-yerr4, E[3]+yerr4, color=line_c.get_color(), alpha =.3)
##Me2
line_c,= ax.plot(nm[4],E[4], color ="cyan", linewidth = 10)
ax.fill_between(nm[4], E[4]-yerr5, E[4]+yerr5, color=line_c.get_color(), alpha =.3)
##Me3
line_c,= ax.plot(nm[5],E[5], color ="orange", linewidth = 10)
ax.fill_between(nm[5], E[5]-yerr[5], E[5]+yerr6, color=line_c.get_color(), alpha =.3)
#ax.fill_between(nm6, E6-.35, E6+.35, color=line_c.get_color(), alpha =.3)

##4 Mes
line_c,= ax.plot(nm[6], E[6], color ="green", linewidth = 10)
ax.fill_between(nm[6], E[6]-yerr7, E[6]+yerr7, color=line_c.get_color(), alpha =.3)

##4 Mes
line_c,= ax.plot(nm[7], E[7], color ="red", linewidth = 10)
ax.fill_between(nm[7], E[7]-yerr8, E[7]+yerr8, color=line_c.get_color(), alpha =.3)

minima = []
for m in E:
    minima.append(m.min()-m.max())


plt.yticks(minima, rotation = 0,fontsize= 48, fontname='Arial', fontstretch='extra-condensed')#, weight='semibold')
plt.xticks(fontsize=48, fontname='Arial', fontstretch='extra-condensed')
#Set color of ticks to match curve
ax.get_yticklabels()[0].set_color("black")
ax.get_yticklabels()[1].set_color("blue")
ax.get_yticklabels()[2].set_color("hotpink")
ax.get_yticklabels()[3].set_color("purple")
ax.get_yticklabels()[4].set_color("cyan")
ax.get_yticklabels()[5].set_color("orange")
ax.get_yticklabels()[6].set_color("green")
ax.get_yticklabels()[7].set_color("red")

ax.tick_params(direction="in", width=6, length =24, pad=15, bottom= True)#labelsize = 56 commented out 
ax.spines["top"].set_linewidth(10)
ax.spines["bottom"].set_linewidth(10)
ax.spines["left"].set_linewidth(10)
ax.spines["right"].set_linewidth(10)
ax.set_xlim(4.5,5.5)
ax.set_ylim(-175,-105)
#ax.set_tick_params(width=5)
#Despine
#sns.despine()


#plt.title("H3 Tail PMF Curves")
#plt.xlabel("Discrete Coordinates (nm)")
#plt.ylabel('Local Minima (kJ $mol^{-1}$)')
plt.savefig('/home/kmpiston/DNA_analy/PMF/2zoom_pmf.png', dpi=300, transparent=True)