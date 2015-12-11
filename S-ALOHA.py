#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: jairo
# @Date:   2015-12-04 22:06:44
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-08 14:31:58
import numpy as np
import argparse
from scipy import misc as sc
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sigma', type=float, required=True,
                        help='Probabilidad de envío')
    parser.add_argument('-t', '--tau', type=float, required=True,
                        help='Probabilidad de reenvío')
    parser.add_argument('-m', '--nodes', type=int, required=True,
                        help='Cantidad de nodos')
    parser.add_argument('-f', '--file', type=argparse.FileType('w'),
                        required=False, default=sys.stdout,
                        help='Archivo de salida (usado con la opción -o)')
    args = parser.parse_args()

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

    # Imprime a la salida (archivo o stdout) con precision de 8 decimales
    np.savetxt(args.file, P, delimiter=',', fmt='%.8f')


if __name__ == '__main__':
    main()
