import matplotlib.pyplot as plt
import numpy as np


up = np.loadtxt("/results_C12/C12mutualist.txt")
vp = np.loadtxt("results_C12/C12prey.txt")
wp = np.loadtxt("results_C12/C12pred.txt")

plt.plot(up)
plt.plot(vp)
plt.plot(wp)
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("p_tot.png")	