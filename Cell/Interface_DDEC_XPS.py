# coding: utf-8

import ase.io
from ase.io import write
from ase import Atoms
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse
import math
import sys
from io import StringIO

def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Create a carbons' XPS spectra from a
    given DDEC6 analysis input file.""")
    # DDEC filename
    parser.add_argument('-i',
                        help="""Input DDEC file for calculating spectra.""")
    # xyz filename
    parser.add_argument('-xyz',
                        help="""Input cell XYZ file of system.""")

    parser.add_argument('-cell',
                        help="""Input cell measurements in Angstroms, separated by comma e.g. 13.4,12.0,2.0""")

    parser.add_argument('-superCell',
                        help="""How many times You want to replicate original Cell for a supercell in each direction for calculating correction,
                         separated by comma e.g. 3,3,1""")

    # Optional args
    parser.add_argument('-c',
                        help="""The index of alkylic carbon in .xyz coordinates (The first atom in .xyz is 0, second 1....)""")
    return parser.parse_args(argv)


""" 
    Function used for calculating binding energies of carbons' electron's for XPS spectra
    
    Arguments:
        Z - 2D array where first array contains charges [e] and second array contains calculated corrections in [e/Angstrom] of atoms
        c - constant, used for shifting alkylic carbon to 285 eV
    Returns:
        Array of binding energies in eVs
"""
def fun(Z, c):
    chg, corr = Z
    return 13.45 * np.asarray(chg) + 14.4 * np.asarray(corr) + c

"""
    Function for estimating the correction electrostatic potential of the other atoms to carbon's BE
    
    Arguments:
        cellDataframe - DataFrame containing the carbon atoms we want to calculate correction 
        cell - ASE Atoms object which contains the supercell used for calculating the correction
 
    Returns:
        1D array of corrections
"""
def CalculateCorrCharge(cellDataframe,cell):
    DDEC = cell.get_initial_charges()
    corrections = []
    for i in cellDataframe.index:
        correction = 0
        distancesFromAtom = cell.get_distances(int(i), np.arange(0,len(cell),1))
        for j in range(len(DDEC)):
            if (i != j):
                correction += float(DDEC[j]) / distancesFromAtom[j]
        corrections += [correction]

    return corrections

"""
    Function for reading in DDEC charge analysis file
    
    Arguments:
        DDECFilename - DDEC charge analysis filename
        
    Returns:
        The lines corresponding to first part of the analysis output as string
"""
def readDDEC(DDECFilename):
    # Read in charges from DDEC
    file = open(args.i, 'r')
    line = file.readline()  # readline
    file.readline()  # Read comment line
    nr_of_atoms = int(line.replace(' ', '').replace('\n', '').replace('\r', ''))  # Remove unnecessary components
    data = str()
    for i in range(nr_of_atoms):  # Read in rows of geometry
        line = file.readline()
        data += line
    return data

"""
    Method for testing giving filename instead of command line excecution
"""
class test:
    def __init__(self, filename, alkylc, xyzFilename, cell, superCell):
        self.i = filename
        self.c = alkylc
        self.xyz = xyzFilename
        self.cell = cell
        self.superCell = superCell

######################### Main excecution starts here ###################################


#args = process_command_line(sys.argv[1:])
args = test("DDEC6_even_tempered_net_atomic_charges.xyz", None, "run_gr-IL_model-20151216-1-13501_mn.xyz", "34.08,34.433,50","3,3,1")

#Parse the cell length
aux = args.cell.split(",")
xlength = float(aux[0]) #A
ylength = float(aux[1]) #A
zlength = float(aux[2]) #A

#Parse the number of cell to replicate
aux = args.superCell.split(",")
cellsInXDirection = int(aux[0])
cellsInYDirection = int(aux[1])
cellsInZDirection = int(aux[2])

# First column is atom type, second x-coord, third y-coord, fourth z-coord, fifth DDEC partial charge
ddecData = readDDEC(args.i)
data1 = StringIO(ddecData)  # String to behave like IO object
coordsCharge = np.loadtxt(data1, usecols=(1, 2, 3, 4))  # Coordinates and charges

#Read in unit cell from xyz and set cells size and charges
unitCell = ase.io.read(args.xyz, format="xyz")
unitCell.set_cell([xlength,ylength,zlength])
unitCell.set_initial_charges(coordsCharge[:, 3])
unitCell.set_tags(np.arange(0,len(unitCell),1)) #Set tags corresponding to atom nr in xyz

#Create supercell from unit cell
superCell = unitCell.repeat((cellsInXDirection,cellsInYDirection,cellsInZDirection))

# Create a dataframe for easier manipulation with unit cell data
df = pd.DataFrame({'Atom': superCell.get_chemical_symbols(), 'x': superCell.get_positions()[:,0],
                   'y': superCell.get_positions()[:,1],
                   'z': superCell.get_positions()[:,2],
                   'tag':superCell.get_tags(),
                   'charge': superCell.get_initial_charges()})

#Set coordinate limits, where should be the middle cell
xMin = math.floor(cellsInXDirection/2)*xlength
xMax = math.ceil(cellsInXDirection/2)*xlength
yMin = math.floor(cellsInYDirection/2)*ylength
yMax = math.ceil(cellsInYDirection/2)*ylength

#Select atoms of middle cell, ordering dataframe by tags
centerCellDf = df[(df['x']>xMin) & (df['x']<xMax) & (df['y']>yMin) & (df['y']<yMax)]
centerCellDf = centerCellDf.sort_values(by='tag')

# Calculate correction to energy from neighbouring atoms' charges to middle cell atoms
centerCellDf['correction']=CalculateCorrCharge(centerCellDf,superCell)


carbons = centerCellDf[centerCellDf['Atom'] == 'C']  # Select carbons from there
carbonIndexes = np.asarray(centerCellDf['tag'])
carbonIndexes = carbonIndexes[centerCellDf['Atom'] == 'C'] #Carbon indexes

if type(args.c) == type(None):
    #If no alcylic carbon given do not fit anything, just calculate BE values using data and auxilliary variable 0
    BEs = fun([list(carbons['charge'])] + [list(carbons['correction'])], 0)

else:
    # Alkylic carbon's energy is 285 eV
    alkyl_carbon = carbons.iloc[int(args.c)]  # Select alkylic carbon
    # Find the constant value
    popt, pcov = curve_fit(fun, [alkyl_carbon['charge'], alkyl_carbon['correction']], 285)
    # Estimate the BE values by using the function and optimized constant
    # print("Constant value is: " + str(popt))
    BEs = fun([list(carbons['charge'])] + [list(carbons['correction'])], *popt)  # Using function to estimate BEs

#Adding binding energies to dataframe and
carbons['BE'] = BEs

# Create arbitrary X-axis
BE_axis = np.linspace(min(BEs) - 1, max(BEs) + 1, 100)  # eV

# Create intenstity Y-axis
arb_intensity = np.zeros(len(BE_axis))

# Way to find indexes of points, where set intensity 1
indexes = np.round((np.sort(BEs) - min(BE_axis)) / (abs(max(BE_axis) - min(BE_axis)) / len(BE_axis)))
for i in indexes:
    arb_intensity[int(i)] += 1

# Plotting and saving data
f1, (ax) = plt.subplots(1, 1, sharey=False, sharex=True, figsize=(8.3 / 2.54, 8.3 / 2.54))
f1.subplots_adjust(hspace=0)

addToFilename = args.i.split(".")[0]

ax.plot(BE_axis, gaussian_filter(arb_intensity, 0.5), 'k')
ax.set_yticklabels([])
ax.tick_params(axis='both', labelsize=9)

ax.set_xlabel('BEs [eV]', fontsize=12)
ax.set_ylabel('Intensity [arb. units]', fontsize=12)

f1.savefig(addToFilename+"_spectra.png", format="png", dpi=300, bbox_inches='tight')
f1.savefig(addToFilename+"_spectra.svg", format="svg")
np.savetxt(addToFilename+"_spectra.csv", list(zip(BE_axis, gaussian_filter(arb_intensity, 2))), delimiter=',')
centerCellDf.to_csv("testing.xyz", columns=["Atom",'x','y','z'], sep="\t")
#print(df)
df = pd.DataFrame(BEs, columns=["BE [eV]"], index=carbonIndexes)
df.index.name = "Atom IX"
df.to_csv(addToFilename+"_BE.csv")