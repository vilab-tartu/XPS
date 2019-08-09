# XPS
Theoretical XPS plotting

DDEC_XPS.y
Create a carbons' XPS spectra from a given DDEC6 analysis input file.
Use example:
python DDEC_XPS.py -i DDEC6_even_tempered_net_atomic_charges.xyz -c 1
Arguments:
-i Input DDEC file for calculating spectra.
-c The index of alkylic carbon in .xyz coordinates (The first atom in .xyz is 0, second 1....)
