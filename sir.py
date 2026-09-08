"""
Introduction to Computational Science
Tijn van Batenburg
2026

SIR model
"""

import numpy as np
from matplotlib import pyplot as plt

def integrate(N, s0, i0, r0, beta, gamma, t0=0, dt=0.01, steps=1000):
    
    s_values = [(s0/N)]
    i_values = [(i0/N)]
    r_values = [(r0/N)]
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

    return np.array(s_values), np.array(i_values), np.array(r_values), np.array(t_values)

def print_plot(s_values, i_values, r_values, t_values, x_points=None, y_points=None):

    plt.plot(t_values, s_values, label="Susceptible")
    plt.plot(t_values, i_values, label="Infected")
    plt.plot(t_values, r_values, label="Recovered")
    plt.xlabel("t")
    plt.ylabel("N")
    plt.title("SIR model")
    plt.legend()
    plt.grid()

    if y_points is not None:
        plt.plot(x_points, y_points, "o")

    plt.show()

# Vaccination strategy 1
def integrate_1(N, s0, i0, r0, beta, gamma, t0=0, dt=0.01, steps=1000):
    
    s_values = [(s0/N)]
    i_values = [(i0/N)]
    r_values = [(r0/N)]
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

        tmp = 30/N * dt if s*N >= 30 else s * dt

        if t_new > 3:
            s_new -= tmp
            r_new += tmp

        s_values.append(s_new)
        i_values.append(i_new)
        r_values.append(r_new)
        t_values.append(t_new)

    return np.array(s_values), np.array(i_values), np.array(r_values), np.array(t_values)





