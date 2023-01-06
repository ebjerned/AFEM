import matplotlib.pyplot as plt
import numpy as np

up = np.loadtxt("results_C12/C12mutualist.txt")
vp = np.loadtxt("results_C12/C12prey.txt")
wp = np.loadtxt("results_C12/C12pred.txt")
t = np.loadtxt("results_C12/C12time.txt")

plt.plot(t,up)
plt.plot(t,vp)
plt.plot(t,wp)
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("results_C12/p_tot.png")

plt.clf()

plt.plot(t,up)
plt.plot(t,vp)
plt.plot(t,wp)
ax = plt.gca()
ax.set_xlim([0, 150])
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("results_C12/p_tot_zoomed.png")
