#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: jairo
# @Date:   2015-12-08 15:14:46
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-09 17:40:12
import argparse
import numpy as np
from scipy import misc as sc
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import colormaps as cmaps


def throughput(stationay_d, M, sigma, tau):
    """Calcula el throughput para el protocolo S-ALOHA.
       stationay_d  --  Distribución estacionaria del modelo
       M  --  Número de nodos en el modelo
       sigma  --  Probabilidad de envío
       tau  --  Probabilidad de reenvío
       """
    P = np.array([0.]*(M+1))
    for k in range(M):
        P[k] = (M-k)*sigma*(1.-sigma)**(M-k-1)*(1.-tau)**k \
               + k*tau*(1.-tau)**(k-1)*(1.-sigma)**(M-k)

    print P
    return np.sum(P*stationay_d)


def gauss(matrix, sol, n):
    """Resuelve una matriz de transición para obtener la distribución
       estacionaria.

       matrix -- La matriz de transición
       sol -- Un vector con las soluciones iniciales
       n -- Número de incógnitas
       """
    out = sol
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
    return out


def SALOHA_gen(M, sigma, tau):
    """Calcula la matriz de transición para la cadena de Markov
       con los parámetros indicados.

       M -- Cantidad de nodos
       sigma -- Probabilidad de envío
       tau -- Probabilidad de reenvío

       return -- Una matriz de (M+1)x(M+1) con las probabilidades de transición
    """
    P = np.array([0.]*(M+1)**2)
    P.shape = (M+1, M+1,)

    tau = tau
    sigma = sigma
    for i in range(M+1):
        for j in range(M+1):
            if j < (i-1):
                pr = 0.
            elif j == (i-1):
                pr = i*tau*(1.-tau)**(i-1)*(1.-sigma)**(M-i)
            elif j == i:
                pr = (1.-(i*tau*(1.-tau)**(i-1)))*(1.-sigma)**(M-i) + \
                    (M-i)*sigma*(1.-sigma)**(M-i-1)*(1.-tau)**i
            elif j == (i+1):
                pr = (M-i)*sigma*(1.-sigma)**(M-i-1)*(1.-(1.-tau)**i)
            elif j > (i+1):
                pr = sc.comb(M-i, j-i)*sigma**(j-i)*(1.-sigma)**(M-j)
            P[i, j] = pr
    return P


def simulate_MC(tr_matrix, steps, initial_st=0):
    """Simula una cadena de Markov dada su matriz de transición.

       tr_matrix -- La matriz de transición
       steps -- Cantidad de transiciones a realizar en la simulación
       initial_st -- Estado inicial, default 0

       return -- Un vector con la distribución estacionaria obtenida al fin de
                 la simulación.
    """
    total_steps = steps
    St = initial_st
    P = np.array([0.]*tr_matrix.shape[0])
    while steps:
        # U está distribuido uniformemente en el intervalo [0, 1)
        U = np.random.random()
        for i in range(tr_matrix.shape[0]):
            # Sumamos las probabilidades desde el elemento 0 hasta i-1
            s = tr_matrix[St][0:i].sum()
            if s < U and U < (s+tr_matrix[St][i]):
                St = i
                P[i] = P[i] + 1
                break
        steps = steps-1

    return P/total_steps


def main():
    np.set_printoptions(precision=7, suppress=True)
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)

    M = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(0.001, 0.9999, num=50)
    Y = np.linspace(0.001, 0.9999, num=50)
    Z = np.array([0.]*X.shape[0]*Y.shape[0])
    Z.shape = (X.shape[0], Y.shape[0])
    print X
    print '--------------------'
    print Y
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            matrix = SALOHA_gen(M, X[i], Y[j])
            # P1 = gauss(matrix, np.array([1./(M+1)]*(M+1)), M+1)
            P2 = simulate_MC(matrix, 100000)
            print 'SALOHA({0},{1})'.format(X[i], Y[j])
            Z[i][j] = throughput(P2, M, X[i], Y[j])

    X, Y = np.meshgrid(X, Y)
    plt.xlabel('sigma')
    plt.ylabel('tau')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmaps.viridis,
                           linewidth=0, antialiased=False, vmin=0.,
                           vmax=0.5)
    ax.set_zlim(0, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
    return

if __name__ == '__main__':
    main()
