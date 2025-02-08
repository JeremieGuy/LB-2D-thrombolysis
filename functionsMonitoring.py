from numpy import *
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pickle
import csv


# initialising the directories to save fluid output
def createRepositoriesFluid(lattice, fluid, clot):

    # Root monitoring directory
    rootDir = "./Monitoring"
    if not os.path.exists(rootDir):
        os.mkdir(rootDir)
        print("Made new root monitoring directory : " + rootDir)

    mainDirTmp = rootDir + "/FF"

    if lattice.branch:
        txt = "_branch=" + str(lattice.branchSize)
    else :
        txt = "_loop"
    mainDirTmp += txt
    
    mainDirTmp += "_viscosity=" + str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    mainDirTmp += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    mainDirTmp += "_it=" + str(lattice.maxIter)

    if not os.path.exists(mainDirTmp):
        os.mkdir(mainDirTmp)
        print("Made new main monitoring directory : " + mainDirTmp)

    # Directory initialisation
    class Directories:
        root = rootDir
        mainDir = mainDirTmp

    return Directories

# initialising the directories to save thrombolysis output
def createRepositoriesThrombolysis(lattice, fluid, clot, tpa, Dir):

    # Root monitoring directory
    rootDir = "./Monitoring"
    if not os.path.exists(rootDir):
        os.mkdir(rootDir)
        print("Made new root monitoring directory : " + rootDir)

    mainDirTmp = rootDir + "/FF"

    if lattice.branch:
        txt = "_branch=" + str(lattice.branchSize)
    else :
        txt = "_loop"
    mainDirTmp += txt
    
    mainDirTmp += "_viscosity=" + str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    mainDirTmp += "_rhoTPA=" + str(tpa.rho_initial)
    mainDirTmp += "_r=" + str(tpa.r)
    mainDirTmp += "_g=" + str(clot.gamma)
    mainDirTmp += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    mainDirTmp += "_it=" + str(lattice.maxIter)

    if not os.path.exists(mainDirTmp):
        os.mkdir(mainDirTmp)
        print("Made new main monitoring directory : " + mainDirTmp)

    # Directory initialisation
    class Directories:
        root = rootDir
        mainDir = mainDirTmp
        clotFrontDir = ""

    # Directory for Clot leftmost value
    if Dir.clotFront:
        Directories.clotFront = Directories.mainDir + "/clotFront"
        if not os.path.exists(Directories.clotFront):
            os.mkdir(Directories.clotFront)
            print("Made new clot leftmost directory : ", Directories.clotFront)

    return Directories

# Drawing a figure of the system, with clot force and acceleration force fields
def plotSystem(main_directory, lattice, bounceback, openPath, clotMask, pulseField):
    
    # Defining the plotting image 
    nx = lattice.nx
    ny = lattice.ny
    tubeSize = lattice.tubeSize
    branchSize = lattice.branchSize

    # open path = 0
    flags_plot = zeros((nx,ny))
    flags_plot[openPath] = 0

    # bounceback = 1
    flags_plot[bounceback] = 1

    # clot = 2
    flags_plot[clotMask] = 2

    # aceleration pulse field = 3
    flags_plot[pulseField] = 3

    # Generating a figure with labels
    flagsname = ["Open Path", "Bounceback","Clot","Acceleration Field"]

    plt.figure(figsize=(7.9,4))

    values = unique(flags_plot.ravel())
    im = plt.imshow(flags_plot.transpose())

    colors = [im.cmap(im.norm(value)) for value in values]
    patches = [mpatches.Patch(color=colors[i], label=flagsname[i] ) for i in range(len(values)) ]

    plt.title("Flags")
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    
    # Saving
    plt.savefig(main_directory + "/system.png",bbox_inches='tight')
    
    # Cleanup
    # plt.show() 
    plt.close()

# Visualising fluid velocity norms
def visualiseFluidVelocity(u):
    plt.clf()
    plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(), cmap=cm.Reds)
    plt.pause(.01)
    plt.cla()

# Visualising tPA density
def visualiseTPADensity(rho):
    plt.clf()
    plt.imshow(rho.transpose(), cmap=cm.Reds)
    plt.pause(.01)
    plt.cla()

# Saving simulation variables to run simulations with an already converged system
def saveVariables(type, lattice, fluid, clot, fin, fout, rho, u):

    # Creating variable storing directory
    varFolder = "./Variables"
    if not os.path.exists(varFolder):
        os.mkdir(varFolder)
        print("Made new variables storing directory : " + varFolder)

    # File for dumping objects containing all the fluid variables for reference
    filename = varFolder + "/" + type + "_"+str(lattice.nx)+"x"+str(lattice.ny)+"_viscosity="
    filename += str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    filename += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    filename += "_it=" + str(lattice.maxIter)

    # Saving variables
    with open(filename + ".pkl", 'wb') as f:
        pickle.dump([fin, fout, rho, u], f)
    
    print("Saved : ", filename)
    
    # Closing the file
    f.close()

# Recovering simulation already converged variables to start the system
def getVariables(type, lattice, fluid, clot, it):
    # Get correct filename
    varFolder = "./Variables"
    filename = varFolder + "/" + type + "_"+str(lattice.nx)+"x"+str(lattice.ny)+"_viscosity=" 
    filename+= str(fluid.viscosity) + "_Rho=" + str(fluid.rho_initial) 
    filename += "_F=" + str(fluid.F_initial) + "_K=" + str(clot.K_initial)
    filename += "_it=" + str(it)

    # Recovering variables
    with open(filename + ".pkl", "rb") as f:  # Python 3: open(..., 'rb')
        fin, fout, rho, u = pickle.load(f)
    print("loaded : ", filename)

    # Closing the file
    f.close()
    
    # Returning variables
    return fin, fout, rho, u
 
# Get clot front coordinate (relative to clot)
def getFrontIndex(K, Clot, clotMask):
    # get clot coordinates
    rows, cols = where(clotMask)
    row_start, row_end = rows.min(), rows.max() + 1
    col_start, col_end = cols.min(), cols.max() + 1

    # get average clot resistance values 
    Kmean = mean(K[0,row_start:row_end, col_start:col_end].transpose(), axis=0)

    # set front threshold
    threshold = 0.1 * Clot.K_initial[0]

    # get corresponding index
    index = argmax(Kmean > threshold)

    return index

# Generating a csv file with desired data
def saveValues(Directory, file, headerX, headerY, values, iterations):
    # Parse data
    rows = zip(iterations, values)
    
    # Generate file in correct directory
    file_name = Directory + file
    with open(file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write the header
        writer.writerow([headerX, headerY])
        
        # Write the data rows
        writer.writerows(rows)

    print(f"Data has been saved to '{file_name}'.")



