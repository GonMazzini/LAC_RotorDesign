# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 16:29:36 2021

@author: Bruger
"""
import numpy as np
# Cl = [1.5, 1.25, 1.35, 0.5] original input 
Cl = [1.3, 1.25, 1.2, 0.5]
# TODO make sure that the AoA correspond to the final cl values
t_c = [24 ,30 ,36 ,48]

import matplotlib.pyplot as plt

# =============================================================================
# plt.scatter(t_c, Cl)
# plt.ylabel('Cl')
# plt.xlabel('t_c')
# plt.grid()
# =============================================================================

p = np.polyfit(t_c, Cl, deg = 5)
poly = np.poly1d(p)

new_x = np.linspace(18,50,40)
new_y = poly(new_x)

plt.plot(t_c, Cl, "o", new_x, new_y)
plt.ylabel('Cl')
plt.xlabel('t_c')
plt.grid()
plt.show()