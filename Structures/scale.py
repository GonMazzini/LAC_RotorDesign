#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:12:32 2020

@author: ozgo
"""
# ============================================================================
import numpy as np
import os
from pathlib import Path


cwd = os.getcwd()
cwd = cwd.replace('\\','/')
path = Path(cwd).parent

# Output file name
f_new_st = Path(str(path) + '\HAWC_inputs\data\DTU_10MW_RWT_Blade_redesign_st.dat')

# ============================================================================
# INPUTS
# ============================================================================
f_original_st = 'st_original.dat'
f_original_c2 = 'c2_original.dat'
f_new_c2 = 'c2_new.dat'
#
# s_f_0 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
# s_f_1 = [10]
# s_f_2 = []
# s_f_4 = []
#
s_f_0 = [8,9,13,14,16]
s_f_1 = [0,2,3,4,5,6,7,17,18]
s_f_2 = [1,15]
s_f_4 = [10,11,12]
# ============================================================================
st_original = np.loadtxt(f_original_st, skiprows=5)
with open(f_original_st) as f:
    a = f.readlines()
col_name = a[3].split()
#
c2_original = np.loadtxt(f_original_c2, usecols=(1,2,3,4,5),comments=';')
curve_original = np.cumsum ( np.linalg.norm (np.diff (c2_original[:,1:4],
                                                      axis=0), axis=1))
#
c2_new = np.loadtxt(f_new_c2, usecols=(1,2,3,4,5),comments=';')
curve_new = np.cumsum ( np.linalg.norm (np.diff (c2_new[:,1:4],
                                                      axis=0), axis=1))
#
s_r = curve_new[-1]/curve_original[-1]
#
st_new = np.zeros_like(st_original)
#
list_0 = [col_name[i] for i in s_f_0]
list_1 = [col_name[i] for i in s_f_1]
list_2 = [col_name[i] for i in s_f_2]
list_4 = [col_name[i] for i in s_f_4]
#
s_r = (97.77-2.8)/(178.3/2-2.8)
#
for i in s_f_0:
    st_new[:,i] = st_original[:,i] * s_r ** 0.0
#
for i in s_f_1:
    st_new[:,i] = st_original[:,i] * s_r ** 1.0
#
for i in s_f_2:
    st_new[:,i] = st_original[:,i] * s_r ** 2.0
#
for i in s_f_4:
    st_new[:,i] = st_original[:,i] * s_r ** 4.0
# ============================================================================
#
b = []
for i in range(3):
    b += a[i]
b += ('    {:7s}    '*19).format(*(col_name))
b += '\n'
b += a[4]
for i, i_t in enumerate(st_new):
    b += ('{:15.6e}').format((st_original[i,0]*s_r))
    b += ('{:15.6e}'*18).format(*(i_t[1:]))
    b+= '\n'
b = ''.join(b)
fp = open(f_new_st,"w")
fp.writelines(b)
fp.close()
