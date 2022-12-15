import matplotlib.pyplot as plt
from dolfin import * 
from fenics import *
mesh = Mesh("circle.xml.gz")
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
element = MixedElement([P1, P1, P1])
#TH = P1 * P1 * P1
W = FunctionSpace(mesh, element)
#W = VectorFunctionSpace(mesh, "Lagrange",1)
t = 0
T = 50
dt = 0.5
delta1 = 1
delta2 = 1
delta3 = 1
alpha = 0.4
beta = 0.8
gamma = 0.8
zeta = 2
L_0 = 0.4
l = 0.6
m =0.12

class InitialConditions(UserExpression):
	def eval(self,values,x):
		values[0] = 0.1 *( 4/15 - 2*10**-7*(x[0]-0.1*x[1]-350)*(x[0]-0.1*x[1]-67) )
		values[1] = 4/15 - 2*10**-7*(x[0]-0.1*x[1]-350)*(x[0]-0.1*x[1]-67)
		values[2] = 22/25 - 3*10**-5*(x[0]-450)-1.2*10**-4*(x[1]-15)
	
	def value_shape(self):
		return (3,)

u = TrialFunction(W)

q = TestFunction(W)


indata = InitialConditions(degree=2)
u0= Function(W)
u0 = interpolate(indata,W)




a1 = u[0]*q[0]*dx + dt*delta1*inner(grad(u[0]),grad(q[0]))*dx 
b1 = u[1]*q[1]*dx + dt*delta2*inner(grad(u[1]),grad(q[1]))*dx
c1 = u[2]*q[2]*dx + dt*delta3*inner(grad(u[2]),grad(q[2]))*dx

#L1 = dt * ( alpha * u0[0] * (1 - u0[0] / (L_0 + l * u0[1])) * q[0] * dx + u0[0] * q[0] * dx)

#L2 = dt * (beta * u0[1] * (1 - u0[1]) - u0[1] * u0[2] / (alpha + u0[1] + m * u0[0])  * q[1] * dx + u0[1] * q[1] * dx)

#L3 = dt * (-gamma * u0[2] + zeta * u0[1] * u0[2] / (alpha + u0[1] + m * u0[0]) * q[2] * dx + u0[2] * q[2] * dx)

L1 = dt*alpha*u0[0]*(1-u0[0]/(L_0+l*u0[1]))*q[0]*dx + u0[0]*q[0]*dx
L2 = dt*( ( beta* u0[1]*(1-u0[1])-u0[2]*u0[1]/(alpha + u0[1] + m * u0[0])))*q[1]*dx + u0[1] * q[1] * dx
L3 = dt*(-gamma*u0[2]+zeta*u0[1]*u0[2]/(alpha+u0[1]+m*u0[0]))*q[2]*dx + u0[2]*q[2]*dx


out_file = File("results/Csolution.pvd")
out_file2 = File("results/test.pvd")

u = Function(W)
u.assign(u0)

a = a1 + b1 + c1
L = L1 + L2 + L3
F = a - L
u_n = Function(W)
u_n.assign(u0)

M0 = u0[0]*dx
M1 = u0[1]*dx
M2 = u0[2]*dx

pop_u = assemble(M0)
pop_v = assemble(M1)
pop_w = assemble(M2)

pu_tot = []
pv_tot = []
pw_tot = []

time = []
while t < T:
	print(t)
	time.append(t)
	u0.assign(u)	
	#A = assemble(a)
	#b = assemble(L)
	#solve(A, u.vector(), b)
	solve(a==L,u)
	#u0.assign(u)
		
	pop_u = assemble(M0)
	pop_v = assemble(M1)
	pop_w = assemble(M2)
	pu_tot.append(pop_u)
	pv_tot.append(pop_v)
	pw_tot.append(pop_w)

	print("u: " + str(pop_u) + ", v: " + str(pop_v) + ", w: " + str(pop_w))
	out_file << (u,t)
	t += dt

out_file2 << u
plt.plot(pu_tot)
plt.plot(pv_tot)
plt.plot(pw_tot)
plt.savefig("p_tot.png")	
plt.clf()
plt.plot(pv_tot, pw_tot)
plt.savefig("vwphase.png")
#solve(F==0,u_n)

