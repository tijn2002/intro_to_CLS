from sir                import *

from matplotlib         import pyplot as plt
from scipy.optimize     import least_squares

import numpy as np
from scipy import fftpack

def main():

    """
    2.1
    """

    # N = 10000
    # s0 = 9999
    # i0 = 1
    # r0 = 0
    # t0 = 0
    # dt = 0.01
    # steps = 10000
    # mu = 0.02
    # beta = 4
    # gamma = 0.8

    # s21, i21, r21, t21 = integrate_bd(N, s0, i0, r0, mu, beta, gamma, t0, dt, steps)

    """
    2.2
    """

    # N = 10000
    # s0 = 9999
    # i0 = 1
    # r0 = 0
    # t0 = 0
    # dt = 0.01
    # steps = 1000
    # beta = 4
    # gamma = 0.8
    # p = 0.3

    # s22, i22, r22, d22, t22 = integrate_d(N, s0, i0, r0, p, beta, gamma, t0, dt, steps)

    # plt.plot(t22, s22, label="Susceptible")
    # plt.plot(t22, i22, label="Infected")
    # plt.plot(t22, r22, label="Recovered")
    # plt.plot(t22, d22, label="Dead")
    # plt.xlabel("t")
    # plt.ylabel("N")
    # plt.title("SIR model")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # ------------------------------------------------------------
    
    # beta = 8
    # gamma = 2
    # steps = 4000

    # birth_rate = 0.4
    # death_rate = 0.6

    # s22, i22, r22, t22, n22 = integrate_bd_sep(N, s0, i0, r0, birth_rate, death_rate, beta, gamma, t0, dt, steps)

    # # print_plot(s22, i22, r22, t22)

    # plt.plot(t22, n22)
    # plt.xlabel("Time")
    # plt.ylabel("Population N")
    # plt.show()

if __name__ == "__main__":
    main()
