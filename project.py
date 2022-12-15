from dolfin import *

mesh = Mesh("circle.xml.gz")

P1 = FiniteElement("Lagrange", mesh.ufl_cell(),1)
TH = MixedElement([P1, P1, P1])
W = FunctionSpace(mesh, TH)


# ALl source terms f,g,q are zero
t = 0
T = 2	
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
		values[0] = 0 
		values[1] = 4/15 - 2*10**-7*(x[0]-0.1*x[1]-350*(x[0]-0.1*x[1]-67))
		values[2] = 22/25 - 3*10**-5*(x[0]-450)-1.2*10**-4*(x[1]-15)
	
	def value_shape(self):
		return (3,)

u_tot = Function(W)
u, v, w = split(u_tot)

u_ntot = Function(W)
u1,v1,w1 = split(u_ntot)

q1,q2,q3 = TestFunctions(W)


indata = InitialConditions(degree=2)
u0= Function(W)
u0 = interpolate(indata,W)



#a = inner(u, q1)*dx + dt*delta1*inner(grad(u), grad(q2))*dx+\
#    inner(v,q2)*dx +dt*delta2*inner(grad(v), grad(q2))*dx+\
#    inner(w,q3)*dx +dt*delta3*inner(grad(w), grad(q3))*dx



#S = Function(W)
S = Expression("dt*alpha*u0[0]*(1-u0[0]/(L_0+1*u0[1]))", dt = dt, alpha = alpha, u0 = u0, L_0=L_0, degree=2)

#S.vector()[:] =dt*alpha*u0.vector()[:]*(1-u0.vector()[:]/(L_0+l*v0.vector()[:]))
#print(S.vector()[:])
#P = Function(W)
P = Expression("dt*(beta*u0[1]*(1-u0[1]+u0[2]/(alpha+u0[1]+m*u0[0])))",dt =dt, alpha = alpha, u0 = u0, beta = beta, m=m, degree=2)
#P.vector()[:] =dt * (beta * v0.vector()[:] * ( 1 - v0.vector()[:]) + v0.vector()[:] * w0.vector()[:] /( alpha + v0.vector()[:]+ m * u0.vector()[:]))

#Q = Function(W)
Q = Expression("dt*(gamma*u0[2]-zeta*u0[1]*u0[2]/(alpha+u0[1]+m*u0[2]))",dt = dt,alpha=alpha, u0=u0, gamma = gamma, zeta = zeta, m=m, degree=2)
#Q.vector()[:] =dt*(gamma*w0.vector()[:]-zeta*v0.vector()[:]*w0.vector()[:]/(alpha+v0.vector()[:] + m*u0.vector()[:]))

#print(a)
print("\n")

#L =S*q1*dx +dot(u0,q1)*dx + P*q2*dx + v0*q2*dx + Q*q3*dx + w0*q3*dx

k = Constant(dt)
a = u*q1*dx + k*delta1*inner(grad(u),grad(q1))*dx+\
    v*q2*dx + k*delta2*inner(grad(v),grad(q2))*dx+\
    w*q3*dx + k*delta3*inner(grad(w),grad(q3))*dx \
    -(S*q1*dx + u0[0]*q1*dx + P*q2*dx + u0[1]*q2*dx + Q*q3*dx + u0[2]*q3*dx)

#a = u*q1*dx + k*delta1*inner(grad(u),grad(q2))*dx+\
#    v*q2*dx + k*delta2*inner(grad(v),grad(q2))*dx+\
#    w*q3*dx + k*delta3*inner(grad(w),grad(q3))*dx 
#L=  -(S*q1*dx + u0[0]*q1*dx + P*q2*dx + u0[1]*q2*dx + Q*q3*dx + u0[2]*q3*dx)
print(a)
u_t = Function(W)


solve(a==0, u_tot)
print("Loop")
while t<T:
	
	solve(a==0,u1)
	un, vn, wn = split(u1)
	u.assign(un)
	v.assign(vn)
	w.assign(wn)
	#L = inner(S,q1)*dx + inner(un,q1)*dx + inner(P,q2)*dx + inner(vn,q2)*dx + inner(Q,q3)*dx + inner(wn,q3)*dx
	#S.vector()[:] =dt*alpha*un.vector()[:]*(1-un.vector()[:]/(L_0+l*vn.vector()[:]))
	#P.vector()[:] =dt*(beta*vn.vector()[:]*(1-vn.vector()[:])+vn.vector()[:]*wn.vector()[:]/(alpha+vn.vector()[:]+m*un.vector()[:]))
	#Q.vector()[:] =dt*(gamma*wn.vector()[:]-zeta*vn.vector()[:]*wn.vector()[:]/(alpha+vn.vector()[:] + m*un.vector()[:]))
	t += dt	

# Set output file
file = File("test.pvd")
file << u1
#Set initial conditions

