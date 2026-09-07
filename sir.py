"""
Introduction to Computational Science
Tijn van Batenburg
2026

SIR model
"""

import numpy as np
from matplotlib import pyplot as pl

def integrate(s0, i0, r0, beta, gamma, set_s_vals=[], set_i_vals=[], set_r_vals=[]):

    t0 = 0
    dt = 0.01
    steps = 2000
    
    s_values = [s0]
    i_values = [i0]
    r_values = [r0]
    t_values = [t0]

    for step in range(steps):

        s = s_values[-1]
        i = i_values[-1]
        r = r_values[-1]
        t = t_values[-1]

        s_next = -beta * s * i
        i_next = beta * s * i - gamma * i
        r_next = gamma * i

        s_new = s + dt * s_next
        i_new = i + dt * i_next
        r_new = r + dt * r_next
        t_new = t + dt

        s_values.append(s_new)
        i_values.append(i_new)
        r_values.append(r_new)
        t_values.append(t_new)

    return s_values, i_values, r_values, t_values

def integrate_set_i(s0, i0, r0, beta, gamma, set_i_vals=None):

    t0 = 0
    dt = 0.01
    steps = len(set_i_vals)
    
    s_values = [s0]
    i_values = [i0]
    r_values = [r0]
    t_values = [t0]

    for step in range(steps - 1):


        s = s_values[-1]
        i = set_i_vals[step]
        r = r_values[-1]
        t = t_values[-1]


        s_next = -beta * s * i
        # i_next = beta * s * i - gamma * i
        r_next = gamma * i

        s_new = s + dt * s_next
        # i_new = i + dt * i_next
        i_new = set_i_vals[step + 1]
        r_new = r + dt * r_next
        t_new = t + dt

        s_values.append(s_new)
        i_values.append(i_new)
        r_values.append(r_new)
        t_values.append(t_new)

    return s_values, i_values, r_values, t_values




