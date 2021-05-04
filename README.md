# XPS
Theoretical XPS plotting

## DDEC_XPS.py
Creates a carbons' XPS spectra from a given DDEC6 analysis input file.
### Arguments
* -i Input DDEC file for calculating spectra.
* -c The index of alkylic carbon in .xyz coordinates (The first atom in .xyz is 0, second 1....)
### Use example:
```python DDEC_XPS.py -i DDEC6_even_tempered_net_atomic_charges.xyz -c 1```

## CELL
Creates a carbons' XPS spectra from a given DDEC6 analysis input file also allows to create into account effects for periodic system

### Interface_DDEC_XPS.py 

#### Arguments

* -i - Input DDEC file for calculating spectra.
* -xyz - Input cell XYZ file of system.
* -cell' - Input cell measurements in Angstroms, separated by comma e.g. 13.4,12.0,2.0
* -superCell - How many times You want to replicate original Cell for a supercell in each direction for calculating correction, separated by comma e.g. 3,3,1
* -slices - Point out where you want to slice the cell to obtain XPS spectra of different parts of cell, separated by comma e.g. 0,20,30,50
* -c - The index of alkylic carbon in .xyz coordinates (The first atom in .xyz is 0, second 1....)

#### Use example

``` python Interface_DDEC_XPS.py -i DDEC6_-05.xyz -xyz Gr_-05_EMImBF4-4050.xyz -cell 34.08,34.433,50 -superCell 3,3,1 -slices 0,20,30,40```
