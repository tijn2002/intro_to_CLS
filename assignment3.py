from sir                import *

from matplotlib         import pyplot as plt
from scipy.optimize     import least_squares

import numpy as np
from scipy import fftpack

def plot_triangle_boundary():

    x = np.linspace(0, 1, 100)
    plt.plot(x, 1 - x, color="grey", linewidth=1.5)
    plt.plot([0, 1], [0, 0], color="grey", linewidth=1.5)
    plt.plot([0, 0], [0, 1], color="grey", linewidth=1.5)

def print_streamplot_mseir():
    
    S, I = np.meshgrid(np.arange(0, 1.01, 0.01), np.arange(0, 1.01, 0.01))

    d = 0.01
    q = 0.0
    epsilon = 0.2
    beta = 0.5
    gamma = 0.1

    E = beta * S * I / (epsilon + d + q)

    dS = d * (1 - S) - beta * S * I
    dI = epsilon * E - (gamma + d + q) * I

    valid = S + I + E <= 1

    dS = np.ma.masked_where(~valid, dS)
    dI = np.ma.masked_where(~valid, dI)

    plt.streamplot(S, I, dS, dI, density=2, maxlength=10, minlength=0.5, color="r")
    plot_triangle_boundary()

    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)

    plt.xlabel("S")
    plt.ylabel("I")
    plt.title(f"MSEIR model")
    plt.grid()
    plt.show()

    return 0

def main():

    """
    3.1.1
    """
    N = 10000
    m0 = 2000
    e0 = 0
    i0 = 10 
    r0 = 0
    s0 = N - m0 - e0 - i0 - r0
    d = 0.01
    q = 0.0
    delta = 0.1
    epsilon = 0.2
    beta = 0.5
    gamma = 0.1

    # m, s, e, i, r, t = integrate_mseir(N, m0, s0, e0, i0, r0, d, q, delta, epsilon, beta, gamma, dt=0.01, steps=20000)

    # plt.plot(t, m, label="Maternally immune")
    # plt.plot(t, s, label="Susceptible")
    # plt.plot(t, e, label="Exposed")
    # plt.plot(t, i, label="Infected")
    # plt.plot(t, r, label="Recovered")
    # plt.xlabel("t")
    # plt.ylabel("fraction of N")
    # plt.title("MSEIR model")
    # plt.legend()
    # plt.grid()
    # plt.show()

    """
    3.1.2
    """
    # print_streamplot_mseir()

    """
    3.2
    """

    d = 0.02
    q = 0.0

    beta0 = 4
    beta1 = 0.5

    m, s, e, i, r, t = integrate_mseir_seasonal(N, m0, s0, e0, i0, r0, d, q, delta, epsilon, beta0, beta1, gamma,
                                                period=5, dt=0.01, steps=10000)

    plt.plot(t, m, label="Maternally immune")
    plt.plot(t, s, label="Susceptible")
    plt.plot(t, e, label="Exposed")
    plt.plot(t, i, label="Infected")
    plt.plot(t, r, label="Recovered")
    plt.xlabel("t")
    plt.ylabel("fraction of N")
    plt.title("MSEIR model with seasonality")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()