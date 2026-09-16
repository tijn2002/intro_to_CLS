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

def integrate_bd(N, s0, i0, r0, mu, beta, gamma, t0=0, dt=0.01, steps=1000):
    
    s_values = [(s0/N)]
    i_values = [(i0/N)]
    r_values = [(r0/N)]
    t_values = [t0]

    for step in range(steps):

        s = s_values[-1]
        i = i_values[-1]
        r = r_values[-1]
        t = t_values[-1]

        s_next = mu - beta * s * i - mu * s
        i_next = beta * s * i - gamma * i - mu * i
        r_next = gamma * i - mu * r

        s_new = s + dt * s_next
        i_new = i + dt * i_next
        r_new = r + dt * r_next
        t_new = t + dt

        s_values.append(s_new)
        i_values.append(i_new)
        r_values.append(r_new)
        t_values.append(t_new)

    return np.array(s_values), np.array(i_values), np.array(r_values), np.array(t_values)

def integrate_bd_sep(N, s0, i0, r0, b_rate, d_rate, beta, gamma, t0=0, dt=0.01, steps=1000):
    
    s_values = [(s0/N)]
    i_values = [(i0/N)]
    r_values = [(r0/N)]
    t_values = [t0]
    N_values = [N * (s_values[0] + i_values[0] + r_values[0])]

    for step in range(steps):

        s = s_values[-1]
        i = i_values[-1]
        r = r_values[-1]
        t = t_values[-1]

        s_next = b_rate - beta * s * i - d_rate * s
        i_next = beta * s * i - gamma * i - d_rate * i
        r_next = gamma * i - d_rate * r

        s_new = s + dt * s_next
        i_new = i + dt * i_next
        r_new = r + dt * r_next
        t_new = t + dt

        s_values.append(s_new)
        i_values.append(i_new)
        r_values.append(r_new)
        t_values.append(t_new)
        N_values.append(N * (s_new + i_new + r_new))

    return np.array(s_values), np.array(i_values), np.array(r_values), np.array(t_values), np.array(N_values)

# Integrate with specific vaccination strategy
def integrate_1(strategy, N, s0, i0, r0, beta, gamma, t0=0, dt=0.01, steps=1000):
    
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

        if strategy == "1":

            if t_new > 3:

                tmp = 30/N * dt if s*N >= 30 else s * dt

                s_new -= tmp
                r_new += tmp

        if strategy == "2":

            if t < 3 <= t_new:
                
                tmp = 0.5 if s >= 0.5 else s

                s_new -= tmp
                r_new += tmp

        s_values.append(s_new)
        i_values.append(i_new)
        r_values.append(r_new)
        t_values.append(t_new)

    return np.array(s_values), np.array(i_values), np.array(r_values), np.array(t_values)





