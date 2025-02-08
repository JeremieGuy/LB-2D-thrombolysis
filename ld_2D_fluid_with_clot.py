# Author : Guy Jérémie
# Date : 08.02.2025

from numpy import *
from functionsLB import *
from functionsMonitoring import *
import time

####################################### Data Load & Save ###########################################

loadData = False
saveData = True

################################### Flow & Geometry Definition #####################################

# Lattice goemetry definition
class Lattice:
    maxIter = 1000                # Max iterations (dt =1)            
    nx, ny = 260, 200               # Number of lattice nodes (dx = 1)
    tubeSize = 21                   # Diameters of the tubes in the system
    branch = False                  # Determines if the system is a loop or with a colateral branch
    branchSize = 21                 # Sets the width of the colateral branch

# Fluid definition
class Fluid:
    viscosity = 0.01                # Kinematic viscosity
    omega = 1 / (3*viscosity+0.5);  # Relaxation parameter
    rho_initial = 2.5               # Inital density of the system
    F_initial = [0,-0.0001]         # Accelerating force F[nx,ny]
    
# Clot definition
class Clot:
    K_initial = [0.001,0.001]       # Initial resisting force of porous region K[2,nx,ny]
    clotSize = 20                   # Size of the clot lenghtwise in a tube section
    coord = [Lattice.nx//2-clotSize//2, Lattice.nx//2+clotSize//2] # Clot coordinates

####################################### Lattice Constants ###########################################

class D2Q9:
    # Velocity directions vectors, D2Q9
    v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
        [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ]) 
    # Directionnal weights, D2Q9
    w = array([1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36]) 
    # Fluid sound velocity adapted to lattice units
    cs2 = 1/3                                                        

################################## Masks ####################################

# Bounceback nodes mask (Loop = False, Branch = True)
bounceback = generateBouncebackMask(Lattice)

# Open path mask
openPath = invert(bounceback)

# Clot mask for clot in the upper tube
clotMask = generateClotMask(Lattice, Clot)

# Force array resistance for porous region
K = generateK(Lattice, Clot, clotMask)

# Clot remaing values mask
KMask = getKMask(Lattice, K)

# acceleration field for fluid aceleration in the lower left tube section
accField = generateAccFieldMask(Lattice)

# Accelerating force values
F = zeros((2,Lattice.nx, Lattice.ny))
F[0,accField] = Fluid.F_initial[0]
F[1,accField] = Fluid.F_initial[1]

##################### Initialising Output Monitoring Functions #####################

# Generating working directories
Directories = createRepositoriesFluid(Lattice, Fluid, Clot)

# Defining geometry type
if Lattice.branch:
    GeometryType = "branch=" + str(Lattice.branchSize) 
else:
    GeometryType = "loop"

# Display current geometry with clot
plotSystem(Directories.mainDir, Lattice, bounceback, openPath, clotMask, accField)

############################# System Initliaization #################################

# Velocity initialization
vel = zeros((2,Lattice.nx, Lattice.ny))

# Density initialization
rho = full((Lattice.nx, Lattice.ny), Fluid.rho_initial)

# initialization of the populations at equilibrium with the given density & velocity.
fin = equilibrium(rho, vel, Lattice, D2Q9)
fout = equilibrium(rho, vel, Lattice, D2Q9)

# Loading already converged variables for faster execution time
if loadData: fin, fout, _, u = getVariables(GeometryType, Lattice, Fluid, Clot, 65000)

################################# Main time loop ######################################

# Monitoring execution time
start_time = time.time()

# main loop
for execTime in range(Lattice.maxIter):

    # Compute macroscopic variables density and velocity.
    rho, u = macroscopic(fin, Lattice, D2Q9)

    # Compute equilibrium.
    feq = equilibrium(rho, u, Lattice, D2Q9)

    # Fluid BGK collision step for open path
    fout[:,openPath] = fin[:,openPath] - Fluid.omega * (fin[:,openPath] - feq[:,openPath])    
    
    # Bounce-back condition 
    for i in range(9):                                  
        fout[i, bounceback] = fin[8-i, bounceback]

    # Forces (acceleration and clot resistance) application
    fout += addForces(rho, u, F, K, Lattice, D2Q9)
    
    # Streaming step for fluid in every direction i=0:8
    fin[0,:,:] = roll(roll(fout[0,:,:],1,axis=0),1,axis=1)      # i = 0
    fin[1,:,:] = roll(fout[1,:,:],1,axis=0)                     # i = 1
    fin[2,:,:] = roll(roll(fout[2,:,:],1,axis=0),-1,axis=1)     # i = 2
    fin[3,:,:] = roll(fout[3,:,:],1,axis=1)                     # i = 3
    fin[4,:,:] = fout[4,:,:]                                    # i = 4
    fin[5,:,:] = roll(fout[5,:,:],-1,axis=1)                    # i = 5
    fin[6,:,:] = roll(roll(fout[6,:,:],-1,axis=0),1,axis=1)     # i = 6
    fin[7,:,:] = roll(fout[7,:,:],-1,axis=0)                    # i = 7
    fin[8,:,:] = roll(roll(fout[8,:,:],-1,axis=0),-1,axis=1)    # i = 8

    # Visualization of the velocity.
    if (execTime%10==0):
        visualiseFluidVelocity(u)

    # Displaying current progress
    print("iteration : " + str(execTime) + "/" + str(Lattice.maxIter), end="\r")
    

# Final execution time
end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

########################### Converged System Saving ############################# 

# Saving converged system to load directly at next run
if saveData : saveVariables(GeometryType, Lattice, Fluid, Clot, fin, fout, rho, u)
