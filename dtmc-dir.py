#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Jairo Sánchez
# @Date:   2015-12-02 12:03:38
# @Last Modified by:   Jairo Sánchez
# @Last Modified time: 2015-12-02 12:09:19


import argparse
import numpy as np
import scipy as sp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mfile', required=False,
                        help='Archivo con la matriz que modela a la cadena')
    return

if __name__ == '__main__':
    main()