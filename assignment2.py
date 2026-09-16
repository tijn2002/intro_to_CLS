from sir                import *

from matplotlib         import pyplot as plt
from scipy.optimize     import least_squares

import numpy as np

def main():

    """
    2.1
    """

    N = 10000
    s0 = 9999
    i0 = 1
    r0 = 0
    t0 = 0
    dt = 0.01
    steps = 8000
    mu = 0.1
    beta = 3
    gamma = 1.2

    # s21, i21, r21, t21 = integrate_bd(N, s0, i0, r0, mu, beta, gamma, t0, dt, steps)

    # print_plot(s21, i21, r21, t21)


    """
    2.2
    """

    beta = 8
    gamma = 2
    steps = 4000

    birth_rate = 0.4
    death_rate = 0.6

    s22, i22, r22, t22, n22 = integrate_bd_sep(N, s0, i0, r0, birth_rate, death_rate, beta, gamma, t0, dt, steps)

    # print_plot(s22, i22, r22, t22)

    plt.plot(t22, n22)
    plt.xlabel("Time")
    plt.ylabel("Population N")
    plt.show()

if __name__ == "__main__":
    main()
