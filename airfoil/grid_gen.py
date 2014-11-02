#! /usr/bin/env python

import sys

if len(sys.argv) != 3:
  print "usage: ",sys.argv[0]," <x-intervals> <filename>"
  sys.exit(1)

N = int(sys.argv[1])
filename = sys.argv[2]

cells = N**2
nodes = (N+1)**2
edges = 2*N*(N-1)
bedges = 4*N

with open(filename,"w") as f:
  # Write header giving the number of entities present in the mesh
  f.write("%s %s %s %s" % (nodes, cells, edges, bedges))

  # Write out cell coordinates
  f.write("\n".join(["%f %f" % (float(i)/N,float(j)/N) for i in range(N+1) for j in range(N+1)]))

  # Write out cell to node mapping
  f.write("\n".join(["%d %d %d %d" % (i, i+1, i+N+2, i+N+1) for i in range(cells)]))

  # Write out cell to edge mapping

    
