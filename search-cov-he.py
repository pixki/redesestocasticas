#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: jairo
# @Date:   2015-12-10 22:29:54
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-10 23:17:30
import numpy as np
import numpy.random as mtrand
from scipy.stats import rv_continuous
import argparse
import sys


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

    def CoV(self):
        a = np.sqrt(2*self.alpha/self.lambda1 + 2*(1-self.alpha)/self.lambda2 -
                    (self.alpha/self.lambda1 + (1-self.alpha)/self.lambda2)**2)
        return a/self.mean()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cov', type=float, required=True,
                        help='Valor de CoV a buscar.')
    parser.add_argument('-o', '--options', type=int, required=False,
                        help='Cuantos conjuntos de valores a obtener',
                        default=1)
    args = parser.parse_args()

    ps = np.linspace(0.1, 0.99999, num=50)
    lambdas = np.linspace(0.5, 80, num=500)
    for p in ps:
        for lambda1 in lambdas:
            for lambda2 in lambdas:
                he = hyperexp(p, lambda1, lambda2)
                good = np.isclose(he.CoV(), args.cov)
                if good:
                    msg = 'CoV encontrado: {3} HE({0},{1}, {2})'
                    print msg.format(p, lambda1, lambda2, he.CoV())
                    args.options = args.options - 1
                    if args.options == 0:
                        sys.exit()
    print 'No se pudo encontrar todas las alternativas requeridas'


if __name__ == '__main__':
    main()
