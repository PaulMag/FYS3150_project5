import matplotlib.pyplot as plt
import numpy as np
from sys import argv

folder = argv[1]

infile = open("data/%s/info.dat" % folder, "r")

info = infile.readline().split(",")
time, n, dim = float(info[0]), int(info[1]), int(info[2])
#dim = 2 # TODO remove this!

names = infile.readline()[:-2].split(",") # (avoid last comma)
N = len(names)

infile.close()

# If you don't want to plot every step:
#resolution = 0.003
#skip = int( round( resolution / (time / n) ) )
skip = 1

positions = np.zeros((N, n, dim))

for i in range(N):
    infile = open("data/%s/obj%d.dat" % (folder,i), "r")

    for j in range(0, n):
        line = infile.readline().split() # numbers divided by spaces
        for k in range(dim):
            positions[i,j,k] = float(line[k])

    infile.close()

plt.figure()
plt.title(folder.split("__")[0] + \
          ", time = %g years, dt = %g years" % (time, time/n))
plt.axis("equal")
plt.xlabel("x [AU]"); plt.ylabel("y [AU]")
plt.grid('on')

for i in range(N):
    plt.plot( positions[i,::skip,0], positions[i,::skip,1] )
#plt.legend( names, loc="best" ) # No need for label in star cluster.

plt.show()
