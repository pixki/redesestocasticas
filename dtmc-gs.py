#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:38
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-06 00:26:57


import argparse
import numpy as np


def valid(matrix):
    """Realiza validaciones en la cadena introducida"""
    if matrix.shape[0] != matrix.shape[1]:
        return (False, 1, 'Rango incorrecto de la matriz')
    for i in range(matrix.shape[0]):
        if matrix[i].sum() != 1.0:
            return (False, 2, 'Probabilidades de transición incorrectas')
    return (True, 0, 'None')


def gauss(matrix, sol, n):
    out = sol
    print 'Initial Solution: {0}'.format(sol)
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
    print "Convergencia alcanzada en {0} iteraciones".format(iterations)
    return out


def main():
    parser = argparse.ArgumentParser(description='Programa para calcular \
        la distribución estacionaria de cadenas de Markov por el método de \
        Gauss-Seidel.')
    parser.add_argument('-m', '--mfile', required=True,
                        help='Archivo con la matriz que modela a la cadena')
    args = parser.parse_args()
    matrix = np.loadtxt(open(args.mfile, "rb"), delimiter=",")
    print 'Resolviendo la matriz: '
    print matrix
    print matrix.shape
    M = matrix.shape[0]-1
    pi = np.array([1./(M+1)]*matrix.shape[0])
    S = gauss(matrix, pi, M+1)
    print "Distribucion estacionaria encontrada: "
    np.set_printoptions(suppress=True)
    print S
    return

if __name__ == '__main__':
    main()
