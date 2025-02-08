# Author : Guy Jérémie
# Date : 08.02.2025

from numpy import *
from functionsLB2 import *
from functionsMonitoring2 import *
import time

####################################### Data Load & Save ###########################################

loadData = True
saveData = False

################################### Flow & Geometry Definition #####################################

# Lattice goemetry definition
class Lattice:
    maxIter = 100              # Max iterations (dt =1)            
    nx, ny = 260, 200               # Number of lattice nodes (dx = 1)
    tubeSize = 21                   # Diameters of the tubes in the system
    branch = False                  # Determins if the system is a loop or with a colateral branch
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
    gamma = 0.5                     # binding proportion

# tPA definition
class TPA:
    rho_initial = 1                 # tPA added concentration
    r = 0.8                           # tPA reaction proportion

########################## Lattice Constants ###########################################

class D2Q9:
    # Velocity directions vectors, D2Q9
    v = array([ [ 1,  1], [ 1,  0], [ 1, -1], [ 0,  1], [ 0,  0],
        [ 0, -1], [-1,  1], [-1,  0], [-1, -1] ]) 
    # Directionnal weights, D2Q9
    w = array([1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36]) 
    # Fluid sound velocity adapted to lattice units
    cs2 = 1/3                          

class D2Q4:
    # Velocity directions vector for tPA, D2Q4
    v = array([[ 1, 0], [ 0, 1], [ 0, -1], [ -1, 0]])
    # Directionnal weighta for tPA, D2Q4
    w = array([1/4, 1/4, 1/4, 1/4])
    # tPA sound velocity adapted to lattice units
    cs2 = 1/2                                    

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
# Dictionnary to generate directories if needed to save data throughout execution
class DirectoryGen:
    clotFront = True

# Generating working directories
Directories = createRepositoriesThrombolysis(Lattice, Fluid, Clot, TPA, DirectoryGen)

# Defining geometry type
if Lattice.branch:
    GeometryType = "branch=" + str(Lattice.branchSize) 
else:
    GeometryType = "loop"

# Display current geometry with clot
plotSystem(Directories.mainDir, Lattice, bounceback, openPath, clotMask, accField)

# Defining output variables
clotFront = []
iterations = []

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

# tPA density initialization
rhoTPA = zeros((Lattice.nx, Lattice.ny))
rhoTPA[1:Lattice.tubeSize+1, Lattice.ny//2] = TPA.rho_initial

# tPA population initialization
tPAin = equilibriumTPA(rhoTPA, u, Lattice, D2Q4)
tPAout = equilibriumTPA(rhoTPA, u, Lattice, D2Q4)

# tPA binded initialization
tPABind = zeros((4,Lattice.nx, Lattice.ny))

################################# Main time loop ######################################

# Monitoring execution time
start_time = time.time()

# main loop
for execTime in range(Lattice.maxIter):

    # Compute macroscopic variables, density and velocity.
    rho, u = macroscopic(fin, Lattice, D2Q9)            # fluid 
    rhoTPA = macroscopicTPA(tPAin)                      # tPA 

    # injecting tPA constantly
    rhoTPA[1:Lattice.tubeSize+1, Lattice.ny//2] = Fluid.rho_initial

    # Compute equilibrium.
    feq = equilibrium(rho, u, Lattice, D2Q9)           # fluid
    tPAeq = equilibriumTPA(rhoTPA, u, Lattice, D2Q4)   # tPA

    # Fluid BGK collision step for open path
    fout[:,openPath] = fin[:,openPath] - Fluid.omega * (fin[:,openPath] - feq[:,openPath])   

    # tPA BGK collision : only where there is no K
    openPathNoK = where(KMask==False, openPath, False)
    tPAout[:,openPathNoK] = tPAin[:,openPathNoK] - Fluid.omega * (tPAin[:,openPathNoK] - tPAeq[:,openPathNoK])    # tPA
    
    # Bounce-back condition 
    for i in range(9):                                  # fluid
        fout[i, bounceback] = fin[8-i, bounceback]
    for i in range(4):                                  # tPA
        tPAout[i, bounceback] = tPAin[3-i, bounceback]
    for i in range(4):                                  # Partial tPA bounceback on clot nodes
        tPAout[i, KMask] = tPAin[3-i, KMask]

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

    # Streaming step for tPA in every direction i=0:4
    tPAin[0,:,:] = roll(tPAout[0,:,:],1,axis=0)                 # i = 0
    tPAin[1,:,:] = roll(tPAout[1,:,:],1,axis=1)                 # i = 1
    tPAin[2,:,:] = roll(tPAout[2,:,:],-1,axis=1)                # i = 2
    tPAin[3,:,:] = roll(tPAout[3,:,:],-1,axis=0)                # i = 3

    # Bind tPA to clot fribrin
    tPABind, tPAin = bindTPA(Clot, tPAin, tPABind, KMask)

    # Dissolve clot
    K, tPABind = dissolveClot(tPABind, K, TPA)

    # print(KMask.shape)
    KMask = getKMask(Lattice, K)

    # liberate remaining binded tPA for empty sites
    tPABind = liberateTPA(tPABind, KMask)
    
    if(execTime%50==0):
        frontIndex = getFrontIndex(K, Clot, clotMask)
        clotFront.append(frontIndex)
        iterations.append(execTime)

    # Visualization of tPA density
    if (execTime%10==0):
        visualiseTPADensity(rhoTPA)

    # Displaying current progress
    print("iteration : " + str(execTime) + "/" + str(Lattice.maxIter), end="\r")
    

# Final execution time
end_time = time.time()
print("Execution time : " + str(end_time-start_time) + " [s]")

######################## Final Iteration Monitoring ########################## 

saveValues(Directories.clotFront, '/clotFront.csv',
            'it', 'pos', clotFront, iterations)


########################### Converged System Saving ############################# 

# Saving converged system to load directly at next run
if saveData : saveVariables(GeometryType, Lattice, Fluid, Clot, fin, fout, rho, u)