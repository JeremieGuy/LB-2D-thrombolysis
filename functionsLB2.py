from numpy import *
from functionsMonitoring import saveKValues

#################### Main Function Definitions ######################################

# Macroscopic variables
def macroscopic(fin, lattice, d2q9):
    rho = sum(fin, axis=0)
    u = zeros((2, lattice.nx, lattice.ny))
    for i in range(9):
        u[0,:,:] += d2q9.v[i,0] * fin[i,:,:]
        u[1,:,:] += d2q9.v[i,1] * fin[i,:,:]
    u /= rho
    return rho, u

# Fluid equilibrium distribution function
def equilibrium(rho, u, lattice, d2q9): 
    usqr = 3/2 * (u[0]**2 + u[1]**2)
    feq = zeros((9,lattice.nx, lattice.ny))
    for i in range(9):
        cu = 3 * (d2q9.v[i,0]*u[0,:,:] + d2q9.v[i,1]*u[1,:,:])
        feq[i,:,:] = rho*d2q9.w[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq

# Macroscopic variable for tPA (only rho)
def macroscopicTPA(tPAin):
    rhoTPA = sum(tPAin, axis=0)
    return rhoTPA

# tPA Equilibrium distribution function
def equilibriumTPA(rhoTPA, u, lattice, d2q4):
    tPAeq = zeros((4,lattice.nx, lattice.ny))
    for i in range(4):
        vu = d2q4.v[i,0]*u[0,:,:] + d2q4.v[i,1]*u[1,:,:]
        tPAeq[i,:,:] = d2q4.w[i]*rhoTPA*(1 + (1/d2q4.cs2)*vu)
    return tPAeq

# Acceleration force and porous region resistance
def addForces(rho, u, F, K, lattice, d2q9):
    FF = zeros((9,lattice.nx, lattice.ny))
    for i in range(9):
        FF[i,:,:] = rho*(d2q9.v[i,0]*(F[0,:,:] - K[0,:,:]*u[0,:,:]) + d2q9.v[i,1]*(F[1,:,:] - K[1,:,:]*u[1,:,:]))
        FF[i,:,:] = FF[i,:,:] * (d2q9.w[i] / d2q9.cs2)
    return FF

# Mask of clot sites
def getKMask(lattice, K):
    # Initialise mask
    KMask = full((lattice.nx, lattice.ny), False)

    # Check where K is null
    condition = K[0] == 0

    # Update accordingly
    KMask = where(condition, KMask, True)

    return KMask

# Binding tPA to fibrin with a gamma factor
def bindTPA(clot, tPAin, tPABind, KMask):
    # Get binded tPA portion
    for i in range(4):
        tPABind[i,KMask] += clot.gamma*tPAin[i,KMask]

    # Calculate remaining free tPA
    tPAin = tPAin - tPABind

    return tPABind, tPAin

def dissolveClot(tPABind, K, tPA):
    #Compute the total amount of population in tPAin
    sumTPABind = sum(tPABind, axis=0)

    # Compute the dissolution amount (considering isotropic clot)
    dissolutionAmount = abs(sumTPABind*tPA.r*K[0,:,:])

    # Dissolving the clot : 
    # 1. Start by applying dissolution
    K_tmp0 = K[0,:,:] - dissolutionAmount
    K_tmp1 = K[1,:,:] - dissolutionAmount

    # 2. Check if K < 1e-7 -> we consider it to be = 0
    isTooSmall0 = K_tmp0 > 1e-7 
    isTooSmall1 = K_tmp1 > 1e-7

    # 3. Adjust accordingly
    K[0,:,:] = where(isTooSmall0, K_tmp0, 0)
    K[1,:,:] = where(isTooSmall1, K_tmp1, 0)

    # 4. Update tPABind quantities
    for i in range(4):
        tPABind[i,:,:] -= tPABind[i,:,:]*tPA.r

    # Return new tPA population and new clot value
    return K, tPABind

# Removing binded tPA where the clot has been dissolved
def liberateTPA(tPABind, KMask):
    # Create mask of empty K
    emptyK = invert(KMask)

    # set binded tPA to 0 where there is no clot K
    for i in range(4):
        tPABind[i,emptyK] = 0

    return tPABind

################################## Masks Functions ####################################

# Bounceback mask
def generateBouncebackMask(lattice):
    nx = lattice.nx
    ny = lattice.ny
    tubeSize = lattice.tubeSize
    branchSize = lattice.branchSize

    bounceback = full((nx, ny),False)
    bounceback[:,0] = True                              # Top border
    bounceback[:,ny-1] = True                           # Bottom border
    bounceback[0,:] = True                              # Left border
    bounceback[nx-1,:] = True                           # Right border
    bounceback[tubeSize+1:(nx-1-tubeSize),              
            tubeSize+1:(ny-1-tubeSize)] = True          # Center Obstacle
    
    if lattice.branch:
        bounceback[tubeSize+1:(nx-1-tubeSize),
        1+2*tubeSize+1:1+2*tubeSize+1 + branchSize] = False   # Branch
    
    return bounceback

# Acceleration field mask
def generateAccFieldMask(lattice):
    pulseField = full((lattice.nx, lattice.ny),False)
    pulseField[1:lattice.tubeSize+1,(lattice.ny//2-10):(lattice.ny//2+11)] = True
    return pulseField

# Initial Clot mask
def generateClotMask(lattice, clot):
    clotMask = full((lattice.nx, lattice.ny),False)
    clotMask[clot.coord[0]:clot.coord[1]+1,1:1+lattice.tubeSize] = True
    return clotMask

# Clot mask through lysis
def generateK(lattice, clot, clotMask):
    K = zeros((2,lattice.nx, lattice.ny))
    K[0,clotMask] = clot.K_initial[0]
    K[1,clotMask] = clot.K_initial[1]
    return K




