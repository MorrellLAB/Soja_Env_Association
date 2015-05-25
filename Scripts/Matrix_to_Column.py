#!/usr/bin/env python
"""Convert the LD matrix output from LDheatmap from the upper triangular matrix
format to a two-column format for easy plotting and calculation of an
exponential regression. Written by Thomas JY Kono"""

import sys
import math

distances = sys.argv[1]
ld_matrix = sys.argv[2]

#   Read in and store the physical distances
phys_dist = []
with open(distances, 'r') as f:
    for line in f:
            phys_dist.append(int(line.strip()))

ldmat = []
#   And do the same for the LD measures
with open(ld_matrix, 'r') as f:
    for line in f:
        ldmat.append(line.strip().split())

#   The LD matrix is a square matrix, so we know its domensions
for i in xrange(0, len(phys_dist)-1):
    for j in xrange(i+1, len(phys_dist)):
        #   get the absolute value of the distance between markers
        d = str(math.fabs(phys_dist[j] - phys_dist[i]))
        #   And the LD measure
        ld = ldmat[i][j]
        #   And print it in two columns
        print d + '\t' + ld
