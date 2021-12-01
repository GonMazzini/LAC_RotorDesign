# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 13:52:29 2021

@author: Mathieu Pellé
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%% For varying wsp

f1 = './HAWC_inputs/DTU_10MW_redesign_flexible_hawc2s.pwr'
f2 = './HAWC_inputs/DTU_10MW_final_ASIER_flexible_hawc2s.pwr'

paths =  [f1, f2]
leg_names = ['DTU 10MW', 'Redesign 1']
col_names = ['V', 'P', 'T', 'CP', 'CT', 'Q', 'Flap M', 'Edge M',\
             'pitch', 'speed', 'tip x', 'typ y', 'tip z', 'J_rot', 'J_dt']
ylab = ['V', 'P [kW]', 'T [kN]', '$C_P$ [-]', '$C_T$ [-]', 'Q', 'Flap M', 'Edge M',\
             'r$\theta$ [deg]', '$\omega$ [rpm]', 'tip x', 'typ y', 'tip z', 'J_rot', 'J_dt']
idx_plot = [1, 2, 3, 4, 8, 9]

for i in range(len(idx_plot)):

    plt.figure()
    for j in range(len(paths)):

        data = pd.read_csv(paths[j], sep='\s+', skiprows = 1, names = col_names)
        x = data['V']
        y = data[col_names[idx_plot[i]]]
        plt.plot(x, y, '-', label=leg_names[j])

    plt.xlabel('$V_{\inf}$ [m/s]')
    plt.ylabel(ylab[idx_plot[i]])
    plt.legend()
    plt.grid()

#%% For varying tsr

path = './V2/results/steady/tsr/V2_hs2.pwr'

R = 95.32742340218736

leg_names = ['DTU 10MW', 'Redesign 1']
col_names = ['V', 'P', 'T', 'CP', 'CT', 'Q', 'Flap M', 'Edge M',\
             'pitch', 'speed', 'tip x', 'typ y', 'tip z', 'J_rot', 'J_dt']
ylab = ['V', 'P [kW]', 'T [kN]', '$C_P$ [-]', '$C_T$ [-]', 'Q', 'Flap M', 'Edge M',\
             'r$\theta$ [deg]', '$\omega$ [rpm]', 'tip x', 'typ y', 'tip z', 'J_rot', 'J_dt']
idx_plot = [1, 2, 3, 4, 8, 9]

for i in range(len(idx_plot)):

    plt.figure()

    data = pd.read_csv(path, sep='\s+', skiprows = 1, names = col_names)
    x = R*data['speed']*np.pi/30/data['V']
    y = data[col_names[idx_plot[i]]]
    plt.plot(x, y, '-')

    plt.xlabel(r'$\lambda$ [m/s]')
    plt.ylabel(ylab[idx_plot[i]])
    plt.grid()
