#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:55
# @Last Modified by:   jairo
# @Last Modified time: 2016-01-15 19:01:28
import numpy as np
from scipy import misc as sc
from scipy.stats import expon
import argparse


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
    parser.add_argument('-S', '--servers', type=int, required=True,
                        help='Cantidad de servidores disponibles')
    parser.add_argument('-l', '--lambd', type=float, required=True,
                        help='Tasa de arribo de nuevos clientes')
    parser.add_argument('-m', '--mu', type=float, required=True,
                        help='Tasa de egreso de clientes')
    args = parser.parse_args()

    a = bcc_recursive(args.servers, args.lambd, args.mu)
    print a, np.sum(a)

    b = bcc_gauss(4, 2, 1)
    print b, np.sum(b)

    c = bcc_sim(4, 2, 1, 100000)
    print c, np.sum(c)

if __name__ == '__main__':
    main()
