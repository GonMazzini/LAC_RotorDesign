# -*- coding: utf-8 -*-
"""An exercise to plot the mean loads for a series of channels and comapre to the theory.
"""

"""Este script compara los resultados de HAWC2S (.pwr) con los steady results de HAWC2 """


import matplotlib.pyplot as plt
import numpy as np
from _loads_utils import load_stats, load_hawc2s


hawc2s_path = 'C:/Users/Bruger/Documents/GitHub/LAC_RotorDesign/Control/DTU_10MW_redesign.pwr'  # path to .pwr or .opt file
stats_path = './res_steady/stats_mean.txt'  # path to mean steady stats

dz_tb = 119  # distance from hub center to tower base [m]
dz_yb = 3.37  # distance from hub center to yaw bearing [m]
geneff = 0.93  # generator/gearboox efficienty [%]
Mgrav = 6200*1.05  # yaw-bearing pitch moment due to gravity [kNm]
#2.7*105520-2.687*446036 #
# define the channels and names to plot
channels = {4: 'Pitch angle [deg]',
            10: 'Rotor speed [rad/s]',
            13: 'Thrust [kN]',
            70: 'Generator torque [Nm]', # was 72
            100: 'Electrical power [W]', # was 102
            61: 'Angle of attack @ 2/3 R [deg]', # was 63 ... and so on
            64: 'Cl @ 2/3 R [-]',
            17: 'Tower-base FA [kNm]',
            18: 'Tower-base SS [kNm]',
            20: 'Yaw-bearing pitch [kNm]',
            21: 'Yaw-bearing roll [kNm]',
            25: 'Shaft torsion [kNm]',
            26: 'Out-of-plane BRM [kNm]',
            27: 'In-plane BRM [kNm]'}
i_wind = 15  # wind channel, needed for plotting versus wind speed

# load the HAWC2 data from the stats file
files, idxs, data = load_stats(stats_path)
wind = data[:, idxs == i_wind]

# load the stuff we need from the HAWC2S .pwr file for the operational data comparisons
h2s_u, h2s_pitch, h2s_rotspd, h2s_paero, h2s_thrust, h2s_aerotrq = load_hawc2s(hawc2s_path)

# loop over each channels and plot the steady state with the theory line
for iplot, (ichan, name) in enumerate(channels.items()):

    # PART 1. Theoretical lines. Taken from HAWC2S pwr file!
    if ichan == 4:  # pitch angle
        u_theory = h2s_u
        theory = h2s_pitch  # directly take pitch angle
    elif ichan == 10:  # rotor speed
        u_theory = h2s_u
        theory = h2s_rotspd  # CORRECT ME!!!
    elif ichan == 13:  # thrust
        u_theory = h2s_u
        theory = h2s_thrust/10**3  # CORRECT ME!!!
    elif ichan == 70:  # generator torque
        u_theory = h2s_u
        theory = h2s_aerotrq*geneff   # CORRECT ME!!!
    elif ichan == 100:  # electrical power
        u_theory = h2s_u
        theory = h2s_paero*geneff  # CORRECT ME!!!

    # extract hawc2 wind and channel to plot from the HAWC2 stats
    h2_wind = data[:, idxs == i_wind]  # wind [m/s]
    HAWC2val = data[:, idxs == ichan]

    # hawc2 channels we need for the theoretical calculations
    h2_thrust = data[:, idxs == 13]  # thrust [kN]
    h2_aero_trq = data[:, idxs == 70] / geneff/10**(3)  # aerodynamic torque [kNm]

    # PART 2. Theoretical lines. Equations in lecture, calculated with hawc2 channels.
    if ichan == 17:  # tower-base fore-aft
        u_theory = h2_wind
        theory = h2_thrust * dz_tb - Mgrav  # tower-base FA is from thrust
    elif ichan == 18:  # tower-base side-side
        u_theory = h2_wind
        theory = h2_aero_trq # CORRECT ME!!!
    elif ichan == 20:  # yaw bearing pitch
        u_theory = h2_wind
        theory = h2_thrust*dz_yb-Mgrav  # CORRECT ME!!!
    elif ichan == 21:  # yaw bearing roll
        u_theory = h2_wind
        theory = h2_aero_trq  # CORRECT ME!!!
    elif ichan == 25:  # shaft torsion
        u_theory = h2_wind
        theory = -h2_aero_trq  # CORRECT ME!!!
    elif ichan == 26:  # blade root out-of-plane-plane moment
        u_theory = h2_wind
        theory = h2_thrust/3  # leave me -- no theory for OoP moment
    elif ichan == 27:  # blade root in-plane moment
        u_theory = h2_wind
        theory = h2_aero_trq/3 # CORRECT ME!!!
    # else:  # no theory
    #     theory = np.full_like(u_theory, np.nan)

    # plot the results
    fig = plt.figure(1 + iplot, figsize=(7, 3), clear=True)
    plt.plot(u_theory, theory, '--', c='0.2')  # theoretical line
    plt.plot(wind, HAWC2val, 'o')  # HAWC2 steady results
    plt.grid('on')
    plt.xlabel('Wind speed [m/s]')
    plt.ylabel(name)
    plt.tight_layout()
    plt.legend(['HAWC2S', 'HAWC2 avg'])
