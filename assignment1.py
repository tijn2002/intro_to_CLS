from sir                import *

from matplotlib         import pyplot as plt
from scipy.optimize     import least_squares

import numpy as np

def print_streamplot():

    S, I = np.meshgrid(np.arange(0, 1.01, 0.01), np.arange(0, 1.01, 0.01))

    beta = 1.5
    gamma = 0.3

    dS = -beta * S * I
    dI = beta * S * I - gamma * I

    valid = S + I <= 1

    dS = np.ma.masked_where(~valid, dS)
    dI = np.ma.masked_where(~valid, dI)

    plt.streamplot(S, I, dS, dI, density=2, maxlength=10, minlength=0.5, color="r")

    beta = 1.5
    gamma = 2

    dS = -beta * S * I
    dI = beta * S * I - gamma * I

    valid = S + I <= 1

    dS = np.ma.masked_where(~valid, dS)
    dI = np.ma.masked_where(~valid, dI)

    plt.streamplot(S, I, dS, dI, density=2, maxlength=10, minlength=0.5, color="b")

    plt.xlabel("S")
    plt.ylabel("I")
    plt.title("SIR model")
    plt.legend()
    plt.grid()
    plt.show()

    return 0

def fit():

    result = least_squares(residuals, x0=[1.0, 0.1], bounds=([0,0],[5,5]))

    beta, gamma = result.x
    mse = np.mean(result.fun ** 2)

    return beta, gamma, mse

def residuals(params):

    i_true = np.array([1, 3, 8, 28, 75, 221, 291, 255, 235, 190, 125, 70, 28])
    i_true = i_true / 763

    beta, gamma = params

    _, i_values, _, _ = integrate(
        763, 762, 1, 0,
        beta=beta,
        gamma=gamma,
        dt=0.01,
        steps=1200
    )

    i_values = i_values[::100][:13]

    return i_values - i_true

def main():

    """
    #1.1
    """

    """
    #1.2
    """

    i_points = (np.array([1, 3, 8, 28, 75, 221, 291, 255, 235, 190, 125, 70, 28])) / 763
    t_points = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

    N = 763
    s0 = 762
    i0 = 1
    r0 = 0
    t0 = 0
    dt = 0.01
    steps = 2000

    beta, gamma, mse = fit()

    # print(beta, gamma, mse)

    # s, i, r, t = integrate(N, s0, i0, r0, beta, gamma, t0, dt, steps)

    # print_plot(s, i, r, t, t_points, i_points)

    """
    #1.3
    """

    # Vaccination strategy 1: Gradual vaccination of 30 children per day starting on day 3

    s, i, r, t = integrate_1(N, s0, i0, r0, beta, gamma, t0, dt, steps)

    print_plot(s, i, r, t, t_points, i_points)

    # Vaccination strategy 2: 

    # Vaccination strategy 3:

if __name__ == "__main__":
    main()





