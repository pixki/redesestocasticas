#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:55
# @Last Modified by:   jairo
# @Last Modified time: 2016-01-13 13:12:10
import numpy as np
from scipy import misc as sc
from scipy.stats import expon
import argparse


def bp_class(S, lambd, mu):
    a = lambd/mu
    pi = np.array([1.]*(S+1))
    j = np.arange(S+1)
    d = np.sum(a**j / sc.factorial(j))
    for i in range(S+1):
        pi[i] = (a**i/sc.factorial(i)) / d

    return pi


def gauss(matrix, sol, n):
    out = sol
    # print 'Initial Solution: {0}'.format(sol)
    iterations = 0
    error = 1000
    while error >= 1e-6:
        old_solution = np.array(out)
        for i in range(n):
            t = np.array(out)
            # Se debe remover el termino del selfloop en el numerador
            t[i] = 0.
            # Se multiplica termino a termino y se suman, \pi_x*P_{xi}
            out[i] = np.sum(t*matrix[:, i]) / (1.-matrix[i][i])

        # Normalizamos el vector para cumplir con la ecuacion \Sum \pi_i=1
        out = out/np.sum(out)
        error = np.sum(np.abs(out - old_solution))
        iterations = iterations + 1
    # print "Convergencia alcanzada en {0} iteraciones".format(iterations)
    return out


def bp_gauss(S, lambd, mu):
    """Genera la matriz con tasas de salida para cada estado de la cadena,
       para que sea resuelta por la funcion gauss
    """
    matrix = np.zeros((S+1, S+1))
    res = np.zeros(S+1)
    for i in range(S+1):
        res[i] = lambd + i*mu
        if i == 0:
            matrix[0][1] = mu
        elif i == S:
            matrix[S][S-1] = lambd
            res[S] = i*mu
        else:
            matrix[i][i-1] = lambd
            matrix[i][i+1] = mu*(i+1)

    return gauss(matrix, res, S+1)


def bp_sim(S, lambd, mu, simtime):
    remaining = simtime
    i = 1   # Estado actual
    ts = 0
    time = np.zeros(S+1)
    while remaining > 0:
        if i == 0:
            T1 = expon.rvs(scale=lambd, size=1)
            T2 = np.inf
        elif i == S:
            T1 = np.inf
            T2 = expon.rvs(scale=i*mu, size=1)
        else:
            T1 = expon.rvs(scale=lambd, size=1)
            T2 = expon.rvs(scale=i*mu, size=1)

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
    a = bp_class(3, 2, 1)
    b = bp_gauss(3, 2, 1)
    c = bp_sim(3, 2, 1, 100000)
    print a, np.sum(a)
    print b, np.sum(b)
    print c, np.sum(c)

if __name__ == '__main__':
    main()
