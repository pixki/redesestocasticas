#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:55
# @Last Modified by:   jairo
# @Last Modified time: 2016-01-18 14:08:29
import numpy as np
from scipy import misc as sc
from scipy.stats import expon
import argparse
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import colormaps as cmaps


def bcc_recursive(S, lambd, mu):
    a = lambd/mu
    pi = np.array([1.]*(S+1))
    j = np.arange(S+1)
    d = np.sum(a**j / sc.factorial(j))
    pi[0] = 1. / d
    for j in range(1, S+1):
        pi[j] = pi[0] * (a**j/sc.factorial(j))
        # print "pi[", j, "]=", a**j, "/", sc.factorial(j), "*", pi[0]

    return pi


def bcc_gauss(S, lambd, mu):
    """Genera la matriz con tasas de salida/entrada para cada estado de
       la cadena, para que sea resuelta por la funcion gauss
    """
    pi = np.array([1./(S+1)]*(S+1))
    error = 1000
    iterations = 0
    while error >= 1e-6:
        old_solution = np.array(pi)
        for i in range(S+1):
            if i == 0:
                pi[0] = mu * pi[1] / (1.*lambd)
            elif i == S:
                pi[S] = (lambd*1.*pi[S-1])/(i*mu*1.)
            else:
                d = (1.*lambd*pi[i-1]+(i+1)*mu*1.*pi[i+1])
                pi[i] = d / (lambd*1. + i*mu*1.)

        pi = pi / np.sum(pi)
        error = np.sum(np.abs(pi - old_solution))
        iterations = iterations + 1

    return pi


def bcc_sim(S, lambd, mu, simtime):
    remaining = simtime
    i = 0   # Estado actual
    ts = 0
    time = np.zeros(S+1)
    while remaining > 0:
        if i == 0:
            T1 = expon.rvs(scale=1./lambd, size=1)
            T2 = np.inf
        elif i == S:
            T1 = np.inf
            T2 = expon.rvs(scale=1./(i*mu), size=1)
        else:
            T1 = expon.rvs(scale=1./lambd, size=1)
            T2 = expon.rvs(scale=1./(i*mu), size=1)

        if np.all(T1 < T2):
            ts = T1
            time[i] = time[i] + ts
            i = i+1
        else:
            ts = T2
            time[i] = time[i] + ts
            i = i-1

        remaining = remaining - ts[0]
        progress = (simtime - remaining) / simtime
        # print "{0}% --> {1} remaining".format(progress*100.0, remaining)

    return time/simtime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--method', type=str, required=True,
                        help='Metodo a usar para resolver el sistema',
                        choices=['gauss', 'simulation', 'recursive'])
    parser.add_argument('-l', '--lambd', type=float, required=True,
                        help='Tasa de arribos, mu se calcula dado el cociente')
    args = parser.parse_args()

    np.set_printoptions(precision=7, suppress=True)
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Con un mayor número de muestras se empieza a 'laggear' el visualizador
    X = np.arange(1, 51)    # Numero de servidores (S=[0, 49])
    Y = np.linspace(0.1, 4.9999, num=50)   # a = lambda/mu
    Z = np.array([0.]*X.shape[0]*Y.shape[0])
    Z.shape = (X.shape[0], Y.shape[0])
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            if 'gauss' in args.method:
                P = bcc_gauss(X[i], args.lambd, args.lambd / Y[i])
            elif 'simulation' in args.method:
                P = bcc_sim(X[i], args.lambd, args.lambd/Y[i], 1000)
            elif 'recursive' in args.method:
                P = bcc_recursive(X[i], args.lambd, args.lambd / Y[i])
            print 'P[S]=', P[-1], ' lambda=', args.lambd, ' mu=',
            print args.lambd/Y[i], ', S=', X[i]

            Z[i][j] = P[-1]

    X, Y = np.meshgrid(X, Y)
    plt.xlabel('S')
    plt.ylabel('a=lambda/mu')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmaps.viridis,
                           linewidth=0, antialiased=True, alpha=1.0,
                           shade=False)
    # ax.set_zlim(0, 1.0)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


if __name__ == '__main__':
    main()
