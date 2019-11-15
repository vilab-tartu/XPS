# coding: utf-8

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import argparse
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
        coordinates - 3D array of atoms' coordinates
        DDEC - 1D array of atoms' DDEC charges
 
    Returns:
        1D array of corrections
"""
def CalculateCorrCharge(coordinates, DDEC):
    x_coords = coordinates[0]
    y_coords = coordinates[1]
    z_coords = coordinates[2]

    corrections = []
    for i in range(len(DDEC)):
        correction = 0
        for j in range(len(DDEC)):
            if (i != j):
                correction += float(DDEC[j]) / (((x_coords[i] - x_coords[j]) ** 2 + (y_coords[i] - y_coords[j]) ** 2 + (
                            z_coords[i] - z_coords[j]) ** 2) ** 0.5)
        corrections += [correction]

    return corrections


"""
    Method for testing giving filename instead of command line excecution
"""
class test:
    def __init__(self, filename, alkylc):
        self.i = filename
        self.c = alkylc

######################### Main excecution starts here ###################################


args = process_command_line(sys.argv[1:])
#args = test("DDEC6_even_tempered_net_atomic_charges.xyz", "0")


file = open(args.i, 'r')
line = file.readline()  # readline
file.readline()  # Read comment line
nr_of_atoms = int(line.replace(' ', '').replace('\n', '').replace('\r', ''))  # Remove unnecessary components
data = str()
for i in range(nr_of_atoms):  # Read in rows of geometry
    line = file.readline()
    data += line

# First column is atom type, second x-coord, third y-coord, fourth z-coord, fifth DDEC partial charge
data1 = StringIO(data)  # String to behave like IO object
atomtype = np.loadtxt(data1, delimiter=' ', usecols=(0), dtype='str')  # Atom types
data1 = StringIO(data)
coordsCharge = np.loadtxt(data1, usecols=(1, 2, 3, 4))  # Coordinates and charges

# Calculate correction to energy from neighbouring atoms' charges

energy_correction = CalculateCorrCharge(np.r_[[coordsCharge[:, 0]], [coordsCharge[:, 1]], [coordsCharge[:, 2]]],
                                        coordsCharge[:, 3])

# Create a dataframe for easier manipulation
df = pd.DataFrame({'Atom': atomtype, 'x': coordsCharge[:, 0], 'y': coordsCharge[:, 1], 'z': coordsCharge[:, 2],
                   'charge': coordsCharge[:, 3], 'correction': energy_correction})
carbons = df[df['Atom'] == 'C']  # Select carbons from there
if type(args.c) == type(None):
    popt, pcov = curve_fit(fun, [carbons['charge'].iloc[0], carbons['correction'].iloc[0]],
                           fun([carbons['charge'].iloc[0]] + [carbons['correction'].iloc[0]], 0))
else:
    # Alkylic carbon's energy is 285 eV
    alkyl_carbon = carbons.iloc[int(args.c)]  # Select alkylic carbon
    # Find the constant value
    popt, pcov = curve_fit(fun, [alkyl_carbon['charge'], alkyl_carbon['correction']], 285)
# Estimate the BE values by using the function
# print("Constant value is: " + str(popt))
BEs = fun([list(carbons['charge'])] + [list(carbons['correction'])], *popt)  # Using function to estimate BEs

# Create arbitrary X-axis
BE_axis = np.linspace(min(BEs) - 1, max(BEs) + 1, 100)  # eV

# Create intenstity Y-axis
arb_intensity = np.zeros(len(BE_axis))

# Way to find indexes of points, where set intensity 1
indexes = np.round((np.sort(BEs) - min(BE_axis)) / (abs(max(BE_axis) - min(BE_axis)) / len(BE_axis)))
for i in indexes:
    arb_intensity[int(i)] = 1

# Plotting and saving data
f1, (ax) = plt.subplots(1, 1, sharey=False, sharex=True, figsize=(8.3 / 2.54, 8.3 / 2.54))
f1.subplots_adjust(hspace=0)

ax.plot(BE_axis, gaussian_filter(arb_intensity, 2), 'k')
ax.set_yticklabels([])
ax.tick_params(axis='both', labelsize=9)

ax.set_xlabel('BEs [eV]', fontsize=12)
ax.set_ylabel('Intensity [arb. units]', fontsize=12)

f1.savefig(args.i.split(".")[0] + ".png", format="png", dpi=300, bbox_inches='tight')
f1.savefig(args.i.split(".")[0] + ".svg", format="svg")
np.savetxt(args.i.split(".")[0] + ".csv", list(zip(BE_axis, gaussian_filter(arb_intensity, 2))), delimiter=',')
