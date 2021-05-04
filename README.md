# XPS
Theoretical XPS plotting

## DDEC_XPS.py
Creates a carbons' XPS spectra from a given DDEC6 analysis input file.
### Arguments
* -i Input DDEC file for calculating spectra.
* -c The index of alkylic carbon in .xyz coordinates (The first atom in .xyz is 0, second 1....)
### Use example:
```python DDEC_XPS.py -i DDEC6_even_tempered_net_atomic_charges.xyz -c 1```

