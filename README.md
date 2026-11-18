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

First, generate a truncated ribosome structure that retains the atoms surrounding the ribosome exit tunnel within an approximately 80 Å radius to reduce computational costs

```
python3 Truncated_6WD4.py
```
Second, use `mole.exe` to perform tunnel searching and extract potential tunnel structures, both with and without specifying the exit-port coordinates. The exit-port coordinates correspond to those of conserved nucleotides in uL23 or eL39 near the tunnel exit, identified through sequence alignment.

```
./mole_Tunnel
```

Then execute `Tunnel_Finalization.py` in Visual Studio to identify the actual exit tunnel for each ribosome, and generate 2D radial plots based on the tunnel centerline coordinates and their corresponding radii.

### MD Simulations

#### Polypeptide Elongation

#### PTC Accessibility
