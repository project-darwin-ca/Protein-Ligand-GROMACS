
# GROMACS-Protein-Ligand

This repository contains resources and scripts for conducting protein-ligand simulations using GROMACS. It provides a step-by-step guide to setting up, running, and analyzing molecular dynamics (MD) simulations.

## Table of Contents
- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Preparation](#1-preparation)
  - [2. Simulation](#2-simulation)
  - [3. Analysis](#3-analysis)
- [File Structure](#file-structure)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Introduction
This project provides a comprehensive workflow for performing protein-ligand MD simulations using GROMACS. It includes scripts for system preparation, simulation execution, and data analysis, helping researchers study the dynamic behavior of protein-ligand complexes.

## Prerequisites
Before using this repository, ensure you have the following software installed:
- GROMACS (version 2020.4 or higher)
- Python (version 3.7 or higher)
- Required Python packages (listed in `requirements.txt`)

## Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/DweipayanG/GROMACS-Protein-Ligand.git
   cd GROMACS-Protein-Ligand
   ```
2. Install the required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### 1. Preparation
Prepare the protein and ligand structures for simulation:
- **Protein preparation**: Clean and preprocess the protein structure using tools like PyMOL or Chimera.
- **Ligand preparation**: Optimize the ligand structure using appropriate software (e.g., Gaussian, Avogadro).

Generate the topology files:
```bash
# Replace protein.pdb and ligand.pdb with your files
./Gromacs\ Codes/generate_topology.sh protein.pdb ligand.pdb
```

### 2. Simulation
Run the molecular dynamics simulations using the prepared files:
```bash
# Minimization
gmx grompp -f EM.mdp -c complex.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# Equilibration
gmx grompp -f NVT.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

gmx grompp -f NPT.mdp -c nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# Production MD
gmx grompp -f MD.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```

### 3. Analysis
Analyze the simulation results using the provided scripts:
```bash
# RMSD calculation
./Gromacs\ Codes/calculate_rmsd.sh md.xtc

# Hydrogen bond analysis
./Gromacs\ Codes/calculate_hbonds.sh md.xtc
```

## File Structure
```
GROMACS-Protein-Ligand/
├── EM.mdp
├── MD.mdp
├── NPT.mdp
├── NVT.mdp
├── ions.mdp
├── cgenff_charmm2gmx_py3_nx.pl
├── sort_mol2_bonds.pl
├── Gromacs Codes/
│   ├── generate_topology.sh
│   ├── calculate_rmsd.sh
│   ├── calculate_hbonds.sh
│   └── ...
├── requirements.txt
└── README.md
```

## Contributing
Contributions are welcome! Please fork this repository, create a new branch, and submit a pull request with your changes. Ensure your code adheres to the existing style and includes appropriate tests.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements
This project was inspired by various GROMACS tutorials and resources. Special thanks to the GROMACS community for their continuous support and development.
