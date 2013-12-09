import numpy as np
from sys import argv
from os import system
import matplotlib.pyplot as plt


### Read files: ###

folder = argv[1]

infile = open("data/%s/info.dat" % folder, "r")

info = infile.readline().split(",")
# objectNo, endTime, timeStep, stepN, dimensionality, snapTime, snapN
N, time, dt, n, dim, snapTime, snapN = \
    int(info[0]), float(info[1]), float(info[2]), \
    int(info[3]), int(info[4]), float(info[5]), int(info[6])

names = infile.readline()[:-2].split(",") # (avoid last comma)
infile.close()

positions = np.zeros((N, snapN, dim))

for i in range(N):
    infile = open("data/%s/obj%d.dat" % (folder,i), "r")

    for j in range(0, snapN):
        line = infile.readline().split() # numbers divided by spaces
        for k in range(dim):
            positions[i,j,k] = float(line[k])

    infile.close()

R0  = 20.0 # "cheating"
eps =  1.5


### 2D trajectory plot: ###
if "2d" in argv or "2D" in argv:

    plt.figure()
    plt.title(folder.split("__")[0] + \
        r", time = %g$\ \tau_{\mathrm{crunch}}$, dt = %g years" % (time, dt))
    plt.axis("equal")
    plt.xlabel("x [ly]"); plt.ylabel("y [ly]")
    plt.grid('on')

    for i in range(N):
        plt.plot( positions[i,:,0], positions[i,:,1] )
    #plt.legend( names, loc="best" ) # No need for label in star cluster.


### 3D movement plot: ###
if "3d" in argv or "3D" in argv:

    from mpl_toolkits.mplot3d import Axes3D

    save = False
    if "film" in argv or "movie" in argv:
        save = True
        system("mkdir data/%s_temp" % folder) # tmp folder for movie frames

    plt.ion()

    fig = plt.figure()
    fig.suptitle("$N = %d$, "        % N + \
        r"$R_0 = %g\ \mathrm{ly}$, " % R0 + \
        r"$T   = %g\ \mathrm{\tau_{crunch}}$, " % time + \
        r"$\mathrm{d}t  = %g\ \mathrm{\tau_{crunch}}$, " % dt + \
        r"$\epsilon = %g\ \mathrm{ly}$"   % eps, fontsize="16")
    ax3 = fig.add_subplot(111, projection='3d')

    ax3.set_xlabel('x [ly]')
    ax3.set_ylabel('y [ly]')
    ax3.set_zlabel('z [ly]')

    dots = []
    for i in range(N):
        # Need to fill dots with something:
        dots.append(ax3.plot([0], [0], [0],'o')[0])
        # Last index [0] because plot returns tuple with one element


    for j in range(0, snapN):
        # This moves around for some reason, not much to do about it:
        ax3.set_title('$t = %.2f$' % (snapTime * j))

        for i in range(N):

            """
            Necessity two have two lines because there is no set_data for
            3 dim data in matplotlib.

            http://matplotlib.org/examples/animation/simple_3danim.html
            """
            dots[i].set_data([positions[i,j,0]],[positions[i,j,1]])
            dots[i].set_3d_properties([positions[i,j,2]])

        ax3.set_xlim([-R0*1.5, R0*1.5])
        ax3.set_ylim([-R0*1.5, R0*1.5])
        ax3.set_zlim([-R0*1.5, R0*1.5])

        plt.draw()
        if save:
            plt.savefig("data/%s_temp/img%06d.png" % (folder, j))

    plt.ioff()

    if save:
        system("ffmpeg -f image2 -r 10 -i data/%s_temp/" % folder + \
               "img%06d.png -vcodec mpeg4 -y " + \
               "data/%s.mp4" % folder)
        # -r is framerate, 6 is good if snapTime=0.02
        system("rm -rf data/%s_temp/" % folder) # delete tmp imgs after film making


### Plot density distribution: ###

def distribution(r):
    n0 = N * 0.001
    r0 = R0 * 2 * N**(-1/3.)
    return n0 / ( 1 + (r / r0)**4 )

if "density" in argv or "dens" in argv:

    r_list = np.zeros(N) # radial distance from origo (center of cluster)
    for i in range(N):
        r_list[i] = np.sqrt(positions[i,-1,0]**2
                          + positions[i,-1,1]**2
                          + positions[i,-1,2]**2)

    dr = 1.0 # bin size
    r_max = np.max(r_list)
    r_n = int(r_max / dr) # number of bins
    r_pts = np.linspace(0, (r_n) * dr, r_n+1)

    r_hist = np.zeros(r_n+1)
    for i in range(N):
        r_hist[ int(r_list[i] / dr) ] += 1 # make a manual histogram

    """
    plt.figure()
    plt.title("$N = %d$, "        % N + \
        r"$R_0 = %g\ \mathrm{ly}$, " % R0 + \
        r"$T   = %g\ \mathrm{\tau_{crunch}}$, " % time + \
        r"$\mathrm{d}t  = %g\ \mathrm{\tau_{crunch}}$, " % dt + \
        r"$\epsilon = %g\ \mathrm{ly}$"   % eps, fontsize="22")
    plt.xlabel("radial distance $[\mathrm{ly}]$", fontsize="26")
    plt.ylabel("star count", fontsize="26")
    plt.bar(r_pts, r_hist, width=dr) # histogram of distribution
    
    plt.plot(r_pts, distribution(r_pts), "r")
    plt.legend(["Theoretical distribution", "Actual counts"], loc="best")
    plt.tight_layout()
    """

    def sphere(r):
        return 4./3. * np.pi * r**3

    for i in range(r_n+1):
        # Divide count by exact volume:
        r_hist[i] /= sphere(r_pts[i] + dr) - sphere(r_pts[i])

    plt.figure()
    plt.title("$N = %d$, "        % N + \
        r"$R_0 = %g\ \mathrm{ly}$, " % R0 + \
        r"$T   = %g\ \mathrm{\tau_{crunch}}$, " % time + \
        r"$\mathrm{d}t  = %g\ \mathrm{\tau_{crunch}}$, " % dt + \
        r"$\epsilon = %g\ \mathrm{ly}$"   % eps, fontsize="22")
    plt.xlabel("radial distance $[\mathrm{ly}]$", fontsize="26")
    plt.ylabel("star density $[1 / \mathrm{ly}^3]$", fontsize="26")
    
    size = 3.0 # wich radius relatice to r_max to plot
    plt.bar(r_pts[0:r_n/size], r_hist[0:r_n/size], width=dr) # histogram of density

    plt.plot(r_pts[0:r_n/size], distribution(r_pts)[0:r_n/size], "r")
    plt.legend(["Theoretical density", "Actual density"], loc="best")
    plt.tight_layout()
    
    print "Average r: ", np.average(r_list)
    print "Stdev of r:", np.std(r_list)

if "2d" in argv or "2D" in argv or "density" in argv or "dens" in argv:
    plt.show()

