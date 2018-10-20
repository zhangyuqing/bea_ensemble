#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os


# In[2]:


## Levels of parameters
method_vec = ["lasso", "rf", "nnet", "svm"]
Nbatch_vec = [2, 3, 4, 5]
batch_mean_vec = [0, 1, 2, 3]
batch_var_vec = [1, 2, 4, 8]
Nsize_vec = [30, 40, 50, 60, 70]


# In[3]:


## Number of batches
sim_short_names = ["nbatch_" + str(ii) + "_" + jj for ii in Nbatch_vec for jj in method_vec]
#print(sim_short_names)
#sim_short_names[0].split("_")

for i in range(len(sim_short_names)):
    splt = sim_short_names[i].split("_")
    os.system("qsub -P combat -o SysEval/logs/" + sim_short_names[i] + ".txt "               + "-e SysEval/logs/" + sim_short_names[i]  + "_err.txt "               + "-N " + sim_short_names[i] +" -l h_rt=72:00:00 -cwd -b y -pe omp 8 "               + "/share/pkg/r/3.4.2/install/bin/Rscript 1_simpipe_syseval.R " + splt[2] + " " + splt[1] + " 1.5 2 FALSE 0 " + splt[0])
    


# In[4]:


## Sample size
sim_short_names = ["size_" + str(ii) + "_" + jj for ii in Nsize_vec for jj in method_vec]

for i in range(len(sim_short_names)):
    splt = sim_short_names[i].split("_")
    os.system("qsub -P combat -o SysEval/logs/" + sim_short_names[i] + ".txt "               + "-e SysEval/logs/" + sim_short_names[i]  + "_err.txt "               + "-N " + sim_short_names[i] +" -l h_rt=72:00:00 -cwd -b y -pe omp 8 "               + "/share/pkg/r/3.4.2/install/bin/Rscript 1_simpipe_syseval.R " + splt[2] + " 3 1.5 2 TRUE " + splt[1] + " " + splt[0])
    


# In[5]:


## Strength of signal (degree of batch effect)
sim_short_names = ["sig_" + str(ii) + "_" + str(jj) + "_" + kk for ii in batch_mean_vec for jj in batch_var_vec for kk in method_vec]
#print(sim_short_names[0].split("_"))

for i in range(len(sim_short_names)):
    splt = sim_short_names[i].split('_')
    os.system("qsub -P combat -o SysEval/logs/" + sim_short_names[i] + ".txt "               + "-e SysEval/logs/" + sim_short_names[i]  + "_err.txt "               + "-N " + sim_short_names[i] +" -l h_rt=72:00:00 -cwd -b y -pe omp 8 "               + "/share/pkg/r/3.4.2/install/bin/Rscript 1_simpipe_syseval.R " + splt[3] + " 3 " + splt[1] + " " + splt[2] + " TRUE 60 " + splt[0])


# In[ ]:




