import matplotlib.pyplot as plt
import random
from dolfin import * 
from fenics import *
mesh = Mesh("meshes/sweden.xml.gz")
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
element = MixedElement([P1, P1, P1])
W = FunctionSpace(mesh, element)

t = 0
T = 1200
dt = 0.5
delta1 = 1
delta2 = 1
delta3 = 1
alpha = 0.4
beta = 1
gamma = 0.8
zeta = 2
L_0 = 0.4
l = 0.4
m = 0.12

class InitialConditions(UserExpression):
	def eval(self,values,x):
		if x[1] > 111:
			values[0] = 1/100
			values[1] = 1/100
			values[2] = 1/100
		else:
			values[0] = 5/1000 * random.uniform(0,1)
			values[1] = 1/2*(1-random.uniform(0,1))
			values[2] = 1/4 + 1/2*random.uniform(0,1)
	
	def value_shape(self):
		return (3,)

u = TrialFunction(W)

q = TestFunction(W)


indata = InitialConditions(degree=2)
u0 = Function(W)
u0 = interpolate(indata,W)




a1 = (1 - dt*alpha) * u[0] * q[0] * dx + dt * delta1 * inner(grad(u[0]), grad(q[0])) * dx 
a2 = (1 - dt*beta)  * u[1] * q[1] * dx + dt * delta2 * inner(grad(u[1]), grad(q[1])) * dx
a3 = (1 + dt*gamma) * u[2] * q[2] * dx + dt * delta3 * inner(grad(u[2]), grad(q[2])) * dx

L1 = - dt * alpha * u0[0] * u0[0] / (L_0 + l * u0[1]) * q[0] * dx + u0[0] * q[0] * dx
L2 = dt * ((-beta * u0[1] * u0[1] - u0[2] * u0[1] / (alpha + u0[1] + m * u0[0]))) * q[1] * dx + u0[1] * q[1] * dx
L3 = dt * ( zeta * u0[1] * u0[2] / (alpha + u0[1] + m * u0[0])) * q[2] * dx + u0[2] * q[2] * dx


out_file = File("results_sweden/SwedenSolution.pvd")

u = Function(W)
u.assign(u0)

a = a1 + a2 + a3
L = L1 + L2 + L3

M0 = u0[0]*dx
M1 = u0[1]*dx
M2 = u0[2]*dx

pop_u = assemble(M0)
pop_v = assemble(M1)
pop_w = assemble(M2)

pu_tot = []
pv_tot = []
pw_tot = []

saving_times = range(0,1200,50) #[0, 50, 100, 600, 1200]

time = []
while t <= T:
	time.append(t)
	u_temp = 0.5*(u0 + u)
	u0.assign(u_temp)	

	solve(a==L,u)

	pop_u = assemble(M0)
	pop_v = assemble(M1)
	pop_w = assemble(M2)
	pu_tot.append(pop_u)
	pv_tot.append(pop_v)
	pw_tot.append(pop_w)

	print("Time: " + str(t) + "u: " + str(pop_u) + ", v: " + str(pop_v) + ", w: " + str(pop_w))
	if t in saving_times:
		out_file << (u,t)
		print("Saved snapshot")
		
	t += dt

plt.plot(pu_tot)
plt.plot(pv_tot)
plt.plot(pw_tot)
plt.legend(["Mutualist", "Prey", "Predator"])
plt.title("Population over time")
plt.xlabel("Time")
plt.ylabel("Population")
plt.savefig("results_sweden/p_tot.png")	

plt.clf()
plt.plot(pv_tot, pw_tot, linewidth=0.2)
plt.title("Phase diagram preys and predators")
plt.xlabel("Preys")
plt.ylabel("Predators")
plt.savefig("results_sweden/vwphase.png")

plt.clf()
plt.plot(pv_tot, pu_tot, linewidth=0.2)
plt.title("Phase diagram preys and mutualists")
plt.xlabel("Preys")
plt.ylabel("Mutualists")
plt.savefig("results_sweden/vuphase.png")

plt.clf()
plt.plot(pw_tot, pu_tot, linewidth=0.2)
plt.title("Phase diagram predators and mutualists")
plt.xlabel("Predators")
plt.ylabel("Mutualists")
plt.savefig("results_sweden/wuphase.png")


