#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:38
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-04 16:33:48


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mfile', required=True,
                        help='Archivo con la matriz de transición.')
    parser.add_argument('-s', '--steps', required=True, type=int,
                        help='Pasos a ejecutar en la simulación.')
    parser.add_argument('-i', '--start', required=True, type=int,
                        help='Estdo inicial de la cadena')
    args = parser.parse_args()
    matrix = np.loadtxt(open(args.mfile, "rb"), delimiter=",")
    print matrix
    r, c, m = valid(matrix)
    if not r:
        print "Cadena no valida: " + m

    steps = args.steps
    St = args.start
    P = np.array([0.]*matrix.shape[0])
    while steps:
        # U está distribuido uniformemente en el intervalo [0, 1)
        U = np.random.random()
        for i in range(matrix.shape[0]):
            # Sumamos las probabilidades desde el elemento 0 hasta i-1
            s = matrix[St][0:i].sum()
            if s < U and U < (s+matrix[St][i]):
                St = i
                P[i] = P[i] + 1
                break
        steps = steps-1

    print P/args.steps
    return

if __name__ == '__main__':
    main()
