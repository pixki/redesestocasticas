#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: jairo
# @Date:   2015-12-04 22:06:44
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-04 23:59:23
import numpy as np
import argparse
from scipy import misc as sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sigma', type=float, required=True,
                        help='Probabilidad de envío')
    parser.add_argument('-t', '--tau', type=float, required=True,
                        help='Probabilidad de reenvío')
    parser.add_argument('-m', '--nodes', type=int, required=True,
                        help='Cantidad de nodos')
    parser.add_argument('-o', '--output', type=str, required=True,
                        choices=['console', 'file'],
                        help='Tipo de salida que ofrece este programa')
    parser.add_argument('-f', '--file', type=str, required=False,
                        help='Archivo de salida (usado con la opción -o)')
    args = parser.parse_args()

    if 'file' in args.output and args.file in '':
        print 'Error: No se ha suministrado archivo de salida'
        exit(-1)

    P = np.array([0.]*(args.nodes+1)**2)
    P.shape = (args.nodes+1, args.nodes+1,)

    tau = args.tau
    sigma = args.sigma
    M = args.nodes
    for i in range(args.nodes+1):
        for j in range(args.nodes+1):
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

    if 'console' in args.output:
        print P
        # Se suma por filas, y se obtiene solamente la ultima columna
        b = np.cumsum(P, axis=1)[:, -1]
        print b
    elif 'file' in args.output:
        np.savetxt(args.file, P, delimiter=',')


if __name__ == '__main__':
    main()
