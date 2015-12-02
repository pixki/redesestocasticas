#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: pixki
# @Date:   2015-11-11 12:07:40
# @Last Modified by:   Jairo Sánchez
# @Last Modified time: 2015-12-02 11:51:32

import numpy as np
from scipy.stats import expon, erlang, rv_continuous
import matplotlib.pyplot as plt
import argparse
import sys

import numpy.random as mtrand


class hyperexp(rv_continuous):

    """An HyperExponential Random Variable
    """

    def __init__(self, alpha=0.5, lambda1=1.0, lambda2=1.0):
        self.alpha = alpha
        self.lambda1 = lambda1
        self.lambda2 = lambda2

    def rvs(self, size=1):
        vsample = np.vectorize(self._single_sample)
        return np.fromfunction(vsample, (size,))

    def _single_sample(self, size):
        U1 = mtrand.random()
        if U1 <= self.alpha:
            scale = self.lambda1
        else:
            scale = self.lambda2
        U2 = mtrand.random()
        return -np.log(U2)/scale

    def pdf(self, x):
        a = self.alpha*self.lambda1*np.exp(self.lambda1*-x)
        b = (1-self.alpha)*self.lambda2*np.exp(self.lambda2*-x)
        return a + b

    def mean(self):
        return (self.alpha / self.lambda1) + ((1-self.alpha) / self.lambda2)

    def standard_dev(self):
        a = (self.alpha/(self.lambda1**2)) + ((1-self.alpha)/(self.lambda2**2))
        return np.sqrt(2*a + self.mean()**2)

    def cdf(self, x):
        a = self.alpha*(-np.exp(self.lambda1*-x))
        b = (1-self.alpha)*(-np.exp(self.lambda2*-x))
        return a + b + 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--stages', type=int, required=False,
                        help='Etapas de la distribución')
    parser.add_argument('-l', '--lambdap', type=float, required=True,
                        nargs='+',
                        help='Parámetro lambda de cada distribución')
    parser.add_argument('-r', '--runs', type=int, required=True,
                        help='Ejecuciones a realizar por cada simulación')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Archivo de salida para la grafica')
    parser.add_argument('-d', '--dist', type=str, required=True,
                        choices=['erlang', 'expon', 'hyperexp'],
                        help='Distribución a emplear para la simulación')
    parser.add_argument('--no-graph', required=False,
                        help='Suprime la salida como gráfica',
                        dest='graph', action='store_false')
    parser.add_argument('--graph', required=False,
                        help='Habilita la salida como gráfica (usar con [-o])',
                        dest='graph', action='store_true')
    parser.set_defaults(graph=True)
    args = parser.parse_args()
    msg = 'Distribución {3} con {0} etapas (lambda={1}) en {2} ejecuciones'
    print msg.format(args.stages, args.lambdap, args.runs, args.dist)
    fig, ax = plt.subplots(1, 1)
    if args.dist in 'erlang':
        if args.stages <= 0:
            print 'Error: se necesita un número válido de etapas'
            sys.exit(1)
        lambdap = args.lambdap[0]
        mean, var, skew, kurt = erlang.stats(args.stages, scale=lambdap,
                                             moments='mvsk')
        print "E[X]={0}, var(X)={1}".format(mean, var)
        x = np.linspace(erlang.ppf(0.00001, args.stages, scale=lambdap),
                        erlang.ppf(0.99999, args.stages, scale=lambdap),
                        num=1000)
        rv = erlang(args.stages, scale=lambdap)
        ax.plot(x, rv.pdf(x), 'r-', lw=5, alpha=0.6, label='Erlang PDF')
        # Generate random numbers with this distribution
        r = erlang.rvs(args.stages, scale=lambdap, size=args.runs)
        ax.hist(r, bins=20, normed=True, histtype='stepfilled', alpha=0.2)
        meanexp = np.mean(r)
        varexp = np.var(r)
        print 'Mediaexperimental: {0} MediaAnalitica: {1}'.format(meanexp,
                                                                  mean)
        print 'Sigma2_exp: {0} Sigma2_a: {1}'.format(varexp, var)
        print 'CoV_exp: {0} CoV_a: {1}'.format(np.sqrt(varexp)/meanexp,
                                               np.sqrt(var)/mean)
    elif args.dist in 'expon':
        lambdap = args.lambdap[0]
        mean, var, skew, kurt = expon.stats(scale=lambdap, moments='mvsk')
        print "E[X]={0}, var(X)={1}".format(mean, var)
        x = np.linspace(expon.ppf(0.00001, scale=lambdap),
                        expon.ppf(0.99999, scale=lambdap),
                        num=1000)
        rv = expon(scale=lambdap)
        ax.plot(x, rv.pdf(x), 'r-', lw=5, alpha=0.6, label='Exponential PDF')
        # Generate random numbers with this distribution
        r = expon.rvs(scale=lambdap, size=args.runs)
        ax.hist(r, bins=20, normed=True, histtype='stepfilled', alpha=0.2)
        meanexp = np.mean(r)
        varexp = np.var(r)
        print 'Mediaexperimental: {0} MediaAnalitica: {1}'.format(meanexp,
                                                                  mean)
        print 'Sigma2_exp: {0} Sigma2_a: {1}'.format(varexp, var)
        print 'CoV_exp: {0} CoV_a: {1}'.format(np.sqrt(varexp)/meanexp,
                                               np.sqrt(var)/mean)
    elif args.dist in 'hyperexp':
        print "HyperExponential RV"
        rv = hyperexp(0.1, args.lambdap[0], args.lambdap[1])
        x = np.linspace(0.00000001, 10.99999, num=1000)
        ax.plot(x, rv.pdf(x), 'r-', lw=5, alpha=0.6, label='HyperExp PDF')
        # ax.plot(x, rv.cdf(x), 'b-', lw=2, alpha=0.6, label='HyperExp CDF')
        r = rv.rvs(size=args.runs)
        ax.hist(r, normed=True, bins=100, range=(0, 11),
                histtype='stepfilled', alpha=0.2)
        meanexp = np.mean(r)
        varexp = np.var(r)
        print 'Mediaexperimental: {0} MediaAnalitica: {1}'.format(meanexp,
                                                                  mean)
        print 'Sigma2_exp: {0} Sigma2_a: {1}'.format(varexp, var)
        print 'CoV_exp: {0} CoV_a: {1}'.format(np.sqrt(varexp)/meanexp,
                                               np.sqrt(var)/mean)
    if args.graph:
        plt.show()

if __name__ == '__main__':
    main()
