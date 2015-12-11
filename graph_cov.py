#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: jairo
# @Date:   2015-12-10 23:57:19
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-11 00:13:03
import numpy as np
from scipy.stats import expon, erlang, rv_continuous
import matplotlib.pyplot as plt
from scipy import misc as sc
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import colormaps as cmaps

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

    def CoV(self):
        a = np.sqrt(2*self.alpha/self.lambda1 + 2*(1-self.alpha)/self.lambda2 -
                    (self.alpha/self.lambda1 + (1-self.alpha)/self.lambda2)**2)
        return a/self.mean()


def main():
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(0.1, 20, num=50)
    Y = np.linspace(0.1, 20, num=50)
    Z = np.array([0.]*X.shape[0]*Y.shape[0])
    Z.shape = (X.shape[0], Y.shape[0])
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            he = hyperexp(0.5, X[i], Y[j])
            Z[i][j] = he.CoV()

    X, Y = np.meshgrid(X, Y)
    plt.xlabel('lambda1')
    plt.ylabel('lambda2')
    ax.set_zlabel('CoV')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmaps.viridis,
                           linewidth=0, antialiased=False, vmin=0., vmax=6)
    ax.set_zlim(0, 6)
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
    return

if __name__ == '__main__':
    main()
