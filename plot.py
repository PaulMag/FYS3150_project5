import numpy as np
from sys import argv
from os import system

folder = argv[1]

infile = open("data/%s/info.dat" % folder, "r")

info = infile.readline().split(",")
time, n, dim = float(info[0]), int(info[1]), int(info[2])
# n = time steps
dt = time / n

names = infile.readline()[:-2].split(",") # (avoid last comma)
N = len(names)
# N = number of objects

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


# Plotting:
import matplotlib.pyplot as plt
"""
plt.figure()
plt.title(folder.split("__")[0] + \
    r", time = %g $\tau_{\mathrm{crunch}}$, dt = %g years" % (time, time/n))
plt.axis("equal")
plt.xlabel("x [l.y.]"); plt.ylabel("y [l.y.]")
plt.grid('on')

for i in range(N):
    plt.plot( positions[i,::skip,0], positions[i,::skip,1] )
#plt.legend( names, loc="best" ) # No need for label in star cluster.

plt.show()
"""

# 3D Plotting:
from mpl_toolkits.mplot3d import Axes3D

save = False
if len(argv) >= 3:
    if argv[2] == "film" or argv[2] == "movie":
        save = True
        system("mkdir data/%s_temp" % folder) # tmp folder for movie frames

plt.ion()

fig = plt.figure()
fig.suptitle(folder.split("__")[0] + \
    r", time = %g $\tau_{\mathrm{crunch}}$, dt = %g $\tau_{\mathrm{crunch}}$" \
    % (time, time/n))
ax3 = fig.add_subplot(111, projection='3d')

ax3.set_xlabel('$x$ [l.y.]')
ax3.set_ylabel('$y$ [l.y.]')
ax3.set_zlabel('$z$ [l.y.]')

dots = []
for i in range(N):
    # Need to fill dots with something:
    dots.append(ax3.plot([0], [0], [0],'o')[0])
    # Last index [0] because plot returns tuple with one element

R0 = np.max(positions[:,0,:]) * 1.05

for j in range(0, n):
    # TODO This moves around for some reason:
    #ax3.set_title('$t = %g$' % (dt * j))
    for i in range(N):

        """
        Necessity two have two lines because there is no set_data for
        3 dim data in matplotlib.

        http://matplotlib.org/examples/animation/simple_3danim.html
        """
        dots[i].set_data([positions[i,j,0]],[positions[i,j,1]])
        dots[i].set_3d_properties([positions[i,j,2]])

    ax3.set_xlim([-R0, R0])
    ax3.set_ylim([-R0, R0])
    ax3.set_zlim([-R0, R0])

    plt.draw()
    if save:
        plt.savefig("data/%s_temp/img%06d.png" % (folder, j))

plt.ioff()

if save:
    system("ffmpeg -f image2 -r 5 -i data/%s_temp/" % folder + \
           "img%06d.png -vcodec mpeg4 -y " + \
           "data/%s.mp4" % folder)
    # -r is framerate
    system("rm -rf data/%s_temp/" % folder) # delete tmp imgs after film making

