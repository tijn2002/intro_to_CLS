"""
Introduction to Computational Science
Tijn van Batenburg
2026

SIR model
"""

import numpy as np
from matplotlib import pyplot as pl

def integrate(s0, i0, r0, beta, gamma):

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

#1.5
s_values, i_values, r_values, t_values = integrate(0.999, 0.001, 0.000, 0.8, 0.4)

# pl.plot(t_values, s_values, label="Susceptible")
# pl.plot(t_values, i_values, label="Infected")
# pl.plot(t_values, r_values, label="Recovered")
# pl.xlabel("t")
# pl.ylabel("N")
# pl.title("SIR model")
# pl.legend()
# pl.grid()
# pl.show()

S, I = np.meshgrid(np.arange(0, 1.01, 0.01), np.arange(0, 1.01, 0.01))

beta = 0.7
gamma = 0.4

dS = -beta * S * I
dI = beta * S * I - gamma * I

valid = S + I <= 1

dS = np.ma.masked_where(~valid, dS)
dI = np.ma.masked_where(~valid, dI)

pl.streamplot(S, I, dS, dI, density=2, maxlength=10, minlength=0.5, color="r")

beta = 1.5
gamma = 2

dS = -beta * S * I
dI = beta * S * I - gamma * I

valid = S + I <= 1

dS = np.ma.masked_where(~valid, dS)
dI = np.ma.masked_where(~valid, dI)

pl.streamplot(S, I, dS, dI, density=2, maxlength=10, minlength=0.5, color="b")

pl.xlabel("S")
pl.ylabel("I")
pl.title("SIR model")
pl.legend()
pl.grid()
pl.show()



