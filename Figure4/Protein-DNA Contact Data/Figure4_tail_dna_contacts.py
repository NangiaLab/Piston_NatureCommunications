#!/usr/bin/env python3
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt 
import MDAnalysis.analysis.distances as distance
import pandas as pd
import seaborn as sns
#import MDAnalysis.core.distances
#from numpy import *
#------------------------------------------------------------------------------------#
# Contact Analysis between two selections,
# 0 to 1 (0=never below threshold, 1=always above threshold)
#------------------------------------------------------------------------------------#

##Paths
#
#path = "/home/kmpiston/Nucleosome_whole_sims/CONTROL/gromacs/" #Control
##short_path = "/home/kmpiston/1histone/charmm-gui-7590708488/gromacs/"#Short
##gggg_path = "/home/kmpiston/2histone/charmm-gui-7592610263/gromacs/" #GGGG
#phos_path = "/home/kmpiston/Nucleosome_whole_sims/PHOS/charmm-gui-0771287668/gromacs_phos/" #Phos
#alk_path = "/home/kmpiston/Nucleosome_whole_sims/ALY/gromacs/" #ALK
#mml_path = '/home/kmpiston/Nucleosome_whole_sims/MML/charmm-gui-0770699177/gromacs_mml/' # MML
#dml_path = '/home/kmpiston/Nucleosome_whole_sims/DML/charmm-gui-0981221004/gromacs/'     #DML
#tml_path = "/home/kmpiston/Nucleosome_whole_sims/TML/charmm-gui-0771900990/gromacs_tml/" #TML
#rep_path = '/home/kmpiston/Nucleosome_whole_sims/MULTI/H3K9K27_TML_H4K16_ALK/charmm-gui-1098645207/gromacs/' #Repressive H3K9TML K27TML
#act_path = '/home/kmpiston/Nucleosome_whole_sims/MULTI/ACTIVE_H3K4K36TML_H3K9K27_MML/'
##Universes
##import tpr and trajectory as universe
#md_u = mda.Universe(path+'md.tpr', path+'center_md.xtc')
##short_u = mda.Universe(short_path+'50.tpr', short_path+'center_short.xtc')
##gggg_u = mda.Universe(gggg_path+'50.tpr', gggg_path+'center_ggg.xtc')
#phos_u = mda.Universe(phos_path+'phos.tpr', phos_path+'center_phos.xtc')
#alk_u = mda.Universe(alk_path+'alk.tpr', alk_path+'center_alk.xtc')
#mml_u = mda.Universe(mml_path+'mml.tpr', mml_path+'center_mml.xtc')
#dml_u = mda.Universe(dml_path+'dml.tpr', dml_path+'center_dml.xtc')
#tml_u = mda.Universe(tml_path+'tml.tpr', tml_path+'center_tml.xtc')
#rep_u = mda.Universe(rep_path+'multi.tpr', rep_path+'center_multi.xtc')
#act_u = mda.Universe(act_path+'active.tpr', act_path+'center_active.xtc')
#
##define your max_distance (cutoff, example 5.0)
##max_distance = *.*
#
#max_distance = 5.0
#timestep = 1
#sel_dna = '(nucleic and name P) and around 35 group chain1'# and name CA' 
#
##H3 Tail Indicies
#md_tail    = 'index 1:585 and not backbone'
##short_tail = 'index 1:429 and not backbone'
##gggg_tail  = 'index 1:499 and not backbone'
#phos_tail  = 'index 1:593 and not backbone'
#alk_tail   = 'index 1:589 and not backbone'
#mml_tail   = 'index 1:588 and not backbone'
#dml_tail   = 'index 1:591 and not backbone'
#tml_tail   = 'index 1:594 and not backbone'
#rep_tail   = 'index 1:603 and not backbone'
#act_tail   = 'index 1:609 and not backbone'
#
#
#u_list = [md_u,phos_u,alk_u,mml_u,dml_u,tml_u,rep_u,act_u]
#tail_list =[md_tail, phos_tail,alk_tail,mml_tail,dml_tail,tml_tail,rep_tail,act_tail]
#
#sys = np.transpose(np.array([u_list,tail_list]))
#
#
#sums = []
#avgs = []
#errs = []
#for u in sys:
#    chain1 = u[0].select_atoms(u[1])
#    chain2 = u[0].select_atoms(sel_dna, chain1=chain1)
#
#    n1 = len(chain1.residues)
#    n2 = len(np.unique(chain2.resids))
#    contact_sum = np.zeros((n1, n2))
#
#    for ts in u[0].trajectory:#[24:28:timestep]:
#        ch1 = chain1.center_of_geometry(compound = 'residues')
#        ch2 = chain2.positions
#        ts_dist = distance.distance_array(ch1,ch2)
#        ts_dist[ts_dist < max_distance] = 1
#        ts_dist[ts_dist > max_distance] = 0
#        contact_sum = ts_dist + contact_sum
#    sums.append(contact_sum)
#    avgs.append((((np.array(contact_sum)).sum(axis=1))/(len(u[0].trajectory)-1)))
#    errs.append((((np.array(contact_sum)).std(axis=1))/np.sqrt(len(u[0].trajectory)-1)))
#
##for n in range(10): #Add in ten 0's for Short systems
##    avgs[1]=np.insert(avgs[1],10,0)
#
##Plotting and DF manipulation
# 
#
#values = np.array(avgs)
#columns = ['WT Control', 'Phos(10&28)', 'K4_Acetyl', 'K4Me1', 'K4Me2', 'K4Me3','K(9&27)Me3','K(4&36)Me3K(9&27)Me1']
#resn =[]
#num = 1 
#chain1 = u_list[0].select_atoms(tail_list[0])
#for n in chain1.residues.resnames:
#    name = str(num)+ ' '+ str(n)
#    resn.append(name)
#    num = num+1
#
#df = pd.DataFrame(values, index=columns, columns=resn)
##Formatting annotation
#label_err = np.array([
#        [errs[0]],
#        [errs[1]],
#        [errs[2]],
#        [errs[3]],
#        [errs[4]],
#        [errs[5]],
#        [errs[6]],
#        [errs[7]]])
##label_mean = np.array([
##        [avgs[0].values],
##        [avgs[1].values],
##        [avgs[2].values],
##        [avgs[3].values],
##        [avgs[4].values],
##        [avgs[5].values],
##        [avgs[6].values],
##        [avgs[7].values]])
#   
##Round to 2 decimals 
#label_err=np.reshape((label_err.round(decimals =2)), (8,38))
##label_mean= label_mean.round(decimals =2)   
###Convert to string
##label_err=label_err.astype(str)
##label_mean=label_mean.astype(str)
##format_array = np.full((label_err.shape),' +/- ')
##
##labels = np.char.add(label_mean,format_array)
##labels = np.char.add(labels,label_err)
##LABELS = np.reshape(labels, (8,38))
#data =np.array(df.values)
#
#labels = (np.asarray(["{1:.1f}\n({0})".format(label_err,data) for label_err, data in zip(label_err.flatten(), data.flatten())])).reshape(8,38)
#
#
#
#df.to_pickle('contacts_df.pkl')
#np.save('contacts_lb',labels)
#np.save('contacts_avg', avgs)
#np.save('contacts_err', errs)
#np.save('conatacts_v',values)

labels   = np.load('/home/kmpiston/DNA_analy/CONTACTS/contacts_lb.npy')
avgs = np.load('/home/kmpiston/DNA_analy/CONTACTS/contacts_avg.npy')
errs = np.load('/home/kmpiston/DNA_analy/CONTACTS/contacts_err.npy')
values   = np.load('/home/kmpiston/DNA_analy/CONTACTS/contacts_v.npy')
df = pd.read_pickle('/home/kmpiston/DNA_analy/CONTACTS/contacts_df.pkl')   
df2 = pd.read_pickle('/home/kmpiston/DNA_analy/CONTACTS/contacts_df.pkl')   


#Delta Map

df2
df2 = df2.transpose()
df2

wt ='WT Control'
s = df2[wt]
df3=df2.sub(s,axis=0).transpose()
df3#
#dfT=df.transpose()
#
#
fig_dims= (70,20)
plt.rcParams['font.size'] = 28
fig, ax = plt.subplots(figsize = fig_dims)
sns.heatmap(df3, 
            #Error 
            #annot=labels, 
            #fmt='', 
            annot=True,
            fmt='.1f',
            annot_kws={'size':36}, 
            robust=True,
            center=0,
            cmap='seismic',
            #vmax=1.25,
            cbar=True,
            linewidths=1, 
            linecolor='white', 
            square=True)

plt.xticks(fontsize=54, fontfamily='arial', fontstretch='condensed')
#
##sns.heatmap(df, annot=True, cmap='plasma', linewidths=1, linecolor='white', vmax=1.25)      
plt.savefig("/home/kmpiston/DNA_analy/CONTACTS/delta_contact_map.png", dpi=300, transparent=True)
