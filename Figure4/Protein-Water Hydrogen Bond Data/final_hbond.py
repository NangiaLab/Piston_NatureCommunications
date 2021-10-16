#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""@author: kmpiston
"""

#Hydrogen Bond Analysis of Histone Protein Tails 

import numpy as np
import MDAnalysis as mda
import pandas as pd
import MDAnalysis.analysis.hbonds as mdh
import matplotlib.pyplot as plt 
import seaborn as sns


#PATHS & DEFINITIONS

#Paths

#path = "/home/kmpiston/Nucleosome_whole_sims/CONTROL/gromacs/" #Control
#phos_path = "/home/kmpiston/Nucleosome_whole_sims/PHOS/charmm-gui-0771287668/gromacs_phos/" #Phosphorylation
#alk_path = "/home/kmpiston/Nucleosome_whole_sims/ALY/gromacs/" #Acetylation
#mml_path = '/home/kmpiston/Nucleosome_whole_sims/MML/charmm-gui-0770699177/gromacs_mml/'#Me1
#dml_path = '/home/kmpiston/Nucleosome_whole_sims/DML/charmm-gui-0981221004/gromacs/' #Me2
#tml_path = "/home/kmpiston/Nucleosome_whole_sims/TML/charmm-gui-0771900990/gromacs_tml/" #Me3
#rep_path = '/home/kmpiston/Nucleosome_whole_sims/MULTI/H3K9K27_TML_H4K16_ALK/charmm-gui-1098645207/gromacs/' #Repressive H3K9TML K27TML
#act_path = '/home/kmpiston/Nucleosome_whole_sims/MULTI/ACTIVE_H3K4K36TML_H3K9K27_MML/' #Active Marks (4 Methylations)
#
##Universes
##import tpr and trajectory as universe
#
#md_u = mda.Universe(path+'md.tpr', path+'center_md.xtc')
#phos_u = mda.Universe(phos_path+'phos.tpr', phos_path+'center_phos.xtc')
#alk_u = mda.Universe(alk_path+'alk.tpr', alk_path+'center_alk.xtc')
#dml_u = mda.Universe(dml_path+'dml.tpr', dml_path+'center_dml.xtc')
#mml_u = mda.Universe(mml_path+'mml.tpr', mml_path+'center_mml.xtc')
#tml_u = mda.Universe(tml_path+'tml.tpr', tml_path+'center_tml.xtc')
#rep_u = mda.Universe(rep_path+'multi.tpr', rep_path+'center_multi.xtc')
#act_u = mda.Universe(act_path+'active.tpr', act_path+'center_active.xtc')
#
##Acceptors and donors to add to standard selection due to Post Translational Modifications
#
#d_ptm=['OT', 'NZ'] 
#a_ptm=['NZ', 'OG', 'O1P', 'OT','O2P']
#
##Initialize variables 
#h_sol=[]                                                                #Store output of MDA Hydrogen Bond analysis
#universes = [md_u,phos_u,alk_u,mml_u,dml_u,tml_u,rep_u,act_u]           #Current Systems
#hb_dfs   = []                                                           #Dataframes all Standard MDA values 
#hb_mdfs  = []                                                           #Melted Dataframe with time, donor, and acceptor
#avg_list = []                                                           #Mean HBonds per Residue per timestep, donors and acceptors combined, water filtered out
#err_list = []                                                           #Standard Error of the Mean for each residue 
#
#
## ANALYSIS
#
##Run standard HBOND analysis from MDA package, put results in dataframe and rearrange to needed format for analysis
#for n in range(len(universes)):
#    h_sol.append((mdh.HydrogenBondAnalysis(universes[n], 'resid 0:37', 'resname TIP3', distance=3.0, angle =150.0, donors=d_ptm, acceptors=a_ptm)))
#    h_sol[n].run()
#    h_sol[n].generate_table()
#    hb_dfs.append(pd.DataFrame.from_records(h_sol[n].table))
#    hb_dfs[n].loc[(hb_dfs[n].donor_resid>38),'donor_resid'] = 50 #All Waters assigned to 50
#    hb_dfs[n].loc[(hb_dfs[n].acceptor_resid>38),'acceptor_resid'] = 50
#    hb_mdfs.append(pd.melt(hb_dfs[n], id_vars=['time','donor_resnm','acceptor_resnm'], value_vars=['donor_resid','acceptor_resid'],value_name='resid'))#Melt Dataframe to stack donors / acceptors for overal count, keep resid
#    avg_list.append((hb_mdfs[n].groupby('resid')['time'].apply(lambda x: x.value_counts().mean())).drop(labels=50))
#    err_list.append((hb_mdfs[n].groupby('resid')['time'].apply(lambda x: x.value_counts().sem())).drop(labels=50))
#
##FORMATTING FOR HEATMAP X AXIS AND LABELS
#
##Pulling and Formatting Resicue names for X axis label
#xnames=[]
#nums = list(range(1,39))
#tail = md_u.atoms.select_atoms('resid 0:37 and name CA').atoms.resnames #Pull string of each residue name
#for n in nums:
#    name =str(n) + ' ' + str(tail[n-1])    
#    xnames.append(name)
##Numpy arrays of average hydrogen bonds and standard error of the mean to change to strings for labels
#bonds=[]
#label_err = []
#for b in avg_list:
#    bonds.append(b.values)
#for e in err_list:
#    label_err.append(e.values)
#bonds = np.array(bonds) #save as np array
#label_err= np.array(label_err)   
#label_err=np.reshape((label_err.round(decimals =2)), (8,38))
#data =np.array(np.reshape((bonds.round(decimals=1)), (8,38)))
#labels = (np.asarray(["{1:.1f}\n({0})".format(label_err,data) for label_err, data in zip(label_err.flatten(), data.flatten())])).reshape(8,38)
#systems = ['Control', 'Phos', 'ALK', 'MML', 'DML', 'TML','Repressive','Active']
#
#
#df = pd.DataFrame(bonds, index=systems, columns=xnames) #final df for visualization


#SAVE / LOAD options us code above and save to save key information, Use load options and comment analysis above for quick plot changes

##Save Commands
#df.to_pickle('df.pkl')
#np.save('lb',labels)
#np.save('avg', avg_list)
#np.save('err', err_list)
#np.save('b',bonds)
#
##Load Commands
#
labels   = np.load('/home/kmpiston/DNA_analy/HBOND/lb.npy')
avg_list = np.load('/home/kmpiston/DNA_analy/HBOND/avg.npy')
err_list = np.load('/home/kmpiston/DNA_analy/HBOND/err.npy')
#bonds   = np.load('/home/kmpiston/DNA_analy/HBOND/b.npy')
df = pd.read_pickle('/home/kmpiston/DNA_analy/HBOND/df.pkl')

#df = df.reindex(['WT Control', 'K4_Acetyl','K4Me1','K4Me2','K4Me3','K(9&27)Me3','Phos(10&28)','K(4&36)Me3K(9&27)Me1'])
df = df.reindex(['Control','ALK','MML','DML','TML','Repressive','Phos','Active'])
abbrev_col = ['A1 ','R2 ','T3 ','K4 ','Q5 ','T6 ','A7 ','R8 ','K9 ','S10','T11'
              , 'G12','G13','K14','A15','P16','R17','K18','Q19','L20','A21',
              'T22','K23','A24','A25','R26','K27','S28','A29','P30','A31',
              'T32','G33','G34','V35','K36','K37','P38']

df.columns = abbrev_col
fig_dims= (11,5)
plt.rcParams['font.size'] = 10
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax = plt.subplots(figsize = fig_dims)
sns.heatmap(df, 
            #Error 
            #annot=labels, 
            #fmt='', 
            annot=True,
            fmt='.1f',
            annot_kws={'size':10,'stretch':'ultra-condensed', 'weight':'semibold'}, 
            robust=True,
            #center=0,
            cmap='Blues',
            cbar_kws={'extend':'both','format':'%.1f','shrink':.47, 'aspect':10,'orientation':'vertical', 'drawedges':False},
            vmax=4.5,
            cbar=False,
            linewidths=1.0, 
            linecolor='white', 
            square=True)


#mono = {'family' : 'monospace'}
plt.xticks(rotation = 90,fontsize=14, fontfamily='monospace', fontstretch='condensed', weight='bold')
ax.axes.get_yaxis().set_ticks([])
ax.tick_params(direction='out', left=False,width=4, length =6)
plt.tight_layout()
#columns = ['Control', 'Phos', 'ALK', 'MML', 'DML', 'TML','Repressive','Active']
#fig_dims= (60,20)
#plt.rcParams['font.size'] = 24
#fig, ax = plt.subplots(figsize = fig_dims)
#sns.heatmap(df, 
#            annot=labels, 
#            fmt='', 
#            annot_kws={'size':20}, 
#            robust=True, 
#            cmap='magma', 
#            cbar=True,
#            cbar_kws={"shrink":.25},
#            linewidths=1, 
#            linecolor='white', 
#            square=True)
#plt.xticks(rotation = 45)
#plt.tight_layout()
plt.savefig("/home/kmpiston/DNA_analy/HBOND/NEW_ORDER_hbond_map.png", dpi=300, transparent=True)

