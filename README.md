# Ribosome Exit Tunnel Extraction Workflow

This is the supplementary code for the paper [Prokaryotic-like ribosome exit tunnels in eukaryotic protists: structural origins and functional consequences](#url_to_be_added).

## Dependencies

- Python 3.9
- [LAMMPS version 29 Aug 2024](https://www.lammps.org/)
- [MOLE2](https://webchem.ncbr.muni.cz/Platform/App/Mole)
- [PyMOL 2.6](https://www.pymol.org/)
  
## Usage

This project involves an automated workflow for tunnel extraction and downstream geometric analysis using Python scripts, as well as coarse-grained molecular dynamics (MD) simulations performed using LAMMPS

- E. coli PDB structure 6WD4 is used here as an illustrative example
  
### Tunnel Extraction & Downstream Geometric Analysis

First, generate a truncated ribosome structure that retains the region surrounding the ribosome exit tunnel within an approximately 80 Å radius to reduce computational costs
'''
python3 Truncated_6WD4.py
'''

### MD Simulations

#### Polypeptide Elongation

#### PTC Accessibility
