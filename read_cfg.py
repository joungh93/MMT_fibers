#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:22:35 2020

@author: jlee
"""


import numpy as np
import pandas as pd


# ----- From target1.cat ----- #
df_cat = pd.read_csv('target1.cat', header=None, skiprows=2, sep='\t', nrows=12192,
                     names=('ra','dec','object','rank','type'))
df_tmp = pd.DataFrame(data = {'target': df_cat.index.values+1})
df_Cat = pd.concat([df_cat, df_tmp], axis=1, sort=False)


# ----- From target1.cfg ----- #
f = open('target1.cfg', 'r')
ll = f.readlines()
f.close()

str_cfg = 'fiber\tra\tdec\tplatex\tplatey\ttarget\trank\n'
str_idx = np.arange(len(ll))[np.array([w == str_cfg for w in ll])]

df_cfg1 = pd.read_csv('target1.cfg', header=None, skiprows=str_idx[0]+2, sep='\t', nrows=300,
                      usecols=(0,1,2,5,6), names=('fiber','ra','dec','target','rank'))

df_cfg2 = pd.read_csv('target1.cfg', header=None, skiprows=str_idx[1]+2, sep='\t', nrows=300,
                      usecols=(0,1,2,5,6), names=('fiber','ra','dec','target','rank'))

df_cfg3 = pd.read_csv('target1.cfg', header=None, skiprows=str_idx[2]+2, sep='\t', nrows=300,
                      usecols=(0,1,2,5,6), names=('fiber','ra','dec','target','rank'))


# ----- Matching target ID ----- #
pd.merge(df_Cat, df_cfg1, on='target').to_csv("./match_cfg1.csv") 
pd.merge(df_Cat, df_cfg2, on='target').to_csv("./match_cfg2.csv")
pd.merge(df_Cat, df_cfg3, on='target').to_csv("./match_cfg3.csv")
