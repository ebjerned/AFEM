import matplotlib.pyplot as plt
import numpy as np

up = np.loadtxt("results_C1/C1mutualist.txt")
vp = np.loadtxt("results_C1/C1prey.txt")
wp = np.loadtxt("results_C1/C1pred.txt")
t = np.loadtxt("results_C1/C1time.txt")

plt.plot(t,up)
plt.plot(t,vp)
plt.plot(t,wp)
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("results_C1/p_tot.png")

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
plt.savefig("results_C1/p_tot_zoomed.png")
