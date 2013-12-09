import matplotlib.pyplot as plt

N = [100, 300, 500, 1000]
r_avg = [12.1028, 14.1343, 10.6845, 12.0686]
r_std = [15.684,  10.4058, 14.7524, 10.4954]

R0 = 20.
time  = 5.
dt = .001
eps = .15
plt.figure()
plt.title( \
    r"$R_0 = %g\ \mathrm{ly}$, " % R0 + \
    r"$T   = %g\ \mathrm{\tau_{crunch}}$, " % time + \
    r"$\mathrm{d}t  = %g\ \mathrm{\tau_{crunch}}$, " % dt + \
    r"$\epsilon = %g\ \mathrm{ly}$"   % eps, fontsize="22")
plt.xlabel("N", fontsize="26")
plt.ylabel("r $[\mathrm{ly}]$", fontsize="26")

plt.plot(N, r_avg, N, r_std)
plt.legend(["Average distance", "Standard deviation"], loc="best")
plt.tight_layout()
plt.show()
