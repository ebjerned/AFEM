from fenics import *

T = 5.0            # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size
eps = 0.01         # diffusion coefficient
K = 10.0           # reaction rate

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





# Read mesh from file
mesh = Mesh('circle.xml.gz')

# Define function space for velocity
V = VectorFunctionSpace(mesh, 'P', 2)

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
W = FunctionSpace(mesh, element)

# Define test functions
q1, q2, q3 = TestFunctions(W)

# Define functions for velocity and concentrations
w =Function(V)
U = Function(W)
u_n = Function(W)
w.vector()[:] = 1.0
# Split system functions to access components
u, v, w = split(U)
u_1, v_1, w_1 = split(u_n)

# Define source terms
f_1 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.1,2)<0.05*0.05 ? 0.1 : 0',
                 degree=1)
f_2 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.1 : 0',
                 degree=1)
f_3 = Constant(0)


class InitialConditions(UserExpression):
	def eval(self,values,x):
		values[0] = 0 
		values[1] = 4/15 - 2*10**-7*(x[0]-0.1*x[1]-350*(x[0]-0.1*x[1]-67))
		values[2] = 22/25 - 3*10**-5*(x[0]-450)-1.2*10**-4*(x[1]-15)
	
	def value_shape(self):
		return (3,)



indata = InitialConditions(degree=2)
u0= Function(W)
u0 = interpolate(indata,W)


S = Expression("dt*alpha*u0[0]*(1-u0[0]/(L_0+1*u0[1]))", dt = dt, alpha = alpha, u0 = u0, L_0=L_0, degree=2)
P = Expression("dt*(beta*u0[1]*(1-u0[1]+u0[2]/(alpha+u0[1]+m*u0[0])))",dt =dt, alpha = alpha, u0 = u0, beta = beta, m=m, degree=2)
Q = Expression("dt*(gamma*u0[2]-zeta*u0[1]*u0[2]/(alpha+u0[1]+m*u0[2]))",dt = dt,alpha=alpha, u0=u0, gamma = gamma, zeta = zeta, m=m, degree=2)




# Define expressions used in variational forms
k = Constant(dt)
K = Constant(K)
eps = Constant(eps)

# Define variational problem
#F = ((u_1 - u_n1) / k)*v_1*dx + dot(w, grad(u_1))*v_1*dx \
#  + eps*dot(grad(u_1), grad(v_1))*dx + K*u_1*u_2*v_1*dx  \
#  + ((u_2 - u_n2) / k)*v_2*dx + dot(w, grad(u_2))*v_2*dx \
#  + eps*dot(grad(u_2), grad(v_2))*dx + K*u_1*u_2*v_2*dx  \
#  + ((u_3 - u_n3) / k)*v_3*dx + dot(w, grad(u_3))*v_3*dx \
#  + eps*dot(grad(u_3), grad(v_3))*dx - K*u_1*u_2*v_3*dx + K*u_3*v_3*dx \
#  - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx

k = Constant(dt)
F = u*q1*dx + k*delta1*inner(grad(u),grad(q2))*dx+\
    v*q2*dx + k*delta2*inner(grad(v),grad(q2))*dx+\
    w*q3*dx + k*delta3*inner(grad(w),grad(q3))*dx \
    -(S*q1*dx + u0[0]*q1*dx + P*q2*dx + u0[1]*q2*dx + Q*q3*dx + u0[2]*q3*dx)

# Create time series for reading velocity data
#timeseries_w = TimeSeries('navier_stokes_cylinder/velocity_series')

# Create VTK files for visualization output
vtkfile_u_1 = File('reaction_system/u_1.pvd')
vtkfile_u_2 = File('reaction_system/u_2.pvd')
vtkfile_u_3 = File('reaction_system/u_3.pvd')

# Create progress bar
#progress = Progress('Time-stepping')
#set_log_level(PROGRESS)

# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Read velocity from file
    # timeseries_w.retrieve(w.vector(), t)

    # Solve variational problem for time step
    solve(F == 0, u0)

    # Save solution to file (VTK)
    _u_1, _u_2, _u_3 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    vtkfile_u_3 << (_u_3, t)

    # Update previous solution
    u_n.assign(u)

    # Update progress bar
    #progress.update(t / T)

# Hold plot
interactive()
