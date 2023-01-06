import matplotlib.pyplot as plt
import numpy as np

up = np.loadtxt("results_sweden/Swedenmutualist.txt")
vp = np.loadtxt("results_sweden/Swedenprey.txt")
wp = np.loadtxt("results_sweden/Swedenpred.txt")
t = np.loadtxt("results_sweden/Swedentime.txt")

plt.plot(t,up)
plt.plot(t,vp)
plt.plot(t,wp)
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("results_sweden/p_tot.png")

plt.clf()

plt.plot(t,up)
plt.plot(t,vp)
plt.plot(t,wp)
ax = plt.gca()
ax.set_xlim([0, 200])
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("results_sweden/p_tot_zoomed.png")
