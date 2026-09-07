from sir.py import integrate
from matplotlib import pyplot as pl
import numpy as np

s_values, i_values, r_values, t_values = integrate(0.999, 0.001, 0.000, 0.8, 0.4)

def print_one():

    pl.plot(t_values, s_values, label="Susceptible")
    pl.plot(t_values, i_values, label="Infected")
    pl.plot(t_values, r_values, label="Recovered")
    pl.xlabel("t")
    pl.ylabel("N")
    pl.title("SIR model")
    pl.legend()
    pl.grid()
    pl.show()

def print_streamplot():

    S, I = np.meshgrid(np.arange(0, 1.01, 0.01), np.arange(0, 1.01, 0.01))

    beta = 1.5
    gamma = 0.3

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

print_streamplot()
