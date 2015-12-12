#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:38
# @Last Modified by:   jairo
# @Last Modified time: 2015-12-11 20:27:39

import sys
import argparse
import numpy as np


def directo(tr_matrix, x0=1.):
    """Calcula la distribución estacionaria de una cadena de Markov
       mediante el método directo. La fórmula para calcular cada término es

              x_{k+1} = x_k - \sum_{i=0; i\neq k}^M  x_i*P_{i,k}
                                       P_{(k+1),k}


       tr_matrix - Matriz de transición a resolver
       x0 - Valor inicial para x_0

       returns - El vector con la distribución estacionaria encontrada
    """
    SD = np.zeros((tr_matrix.shape[0],), dtype=np.float64)
    SD[0] = x0
    for i in range(tr_matrix.shape[0] - 1):
        Px = np.array(tr_matrix[:, i])
        Px[i+1] = 0.
        SD[i+1] = (SD[i] - np.sum(Px*SD))/tr_matrix[i+1, i]

    # Se regresa el vector normalizado
    return SD/np.sum(SD)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mfile', required=False,
                        type=argparse.FileType('r'), default=sys.stdin,
                        help='Archivo con la matriz que modela a la cadena')
    args = parser.parse_args()
    matrix = np.loadtxt(args.mfile, delimiter=",")
    np.set_printoptions(suppress=True, precision=6)
    print directo(matrix, 1)
    return

if __name__ == '__main__':
    main()
