
# Code Folder

This folder contains scripts that are essential for the preparation and processing of molecular dynamics simulations using GROMACS.

## Files

1. **cgenff_charmm2gmx_py3_nx1.py**
   - **Description**: This Python script converts CHARMM General Force Field (CGenFF) parameters to a format compatible with GROMACS.
   - **Usage**: 
     ```sh
     python cgenff_charmm2gmx_py3_nx1.py RESNAME drug.mol2 drug.str charmm36.ff
     ```
   - **Inputs**:
     1. `RESNAME`: Residue name found in the RESI entry in the CHARMM stream file.
     2. `drug.mol2`: The .mol2 file you supplied to ParamChem.
     3. `drug.str`: The CHARMM stream file with topology and parameters generated by ParamChem.
     4. `charmm36.ff`: The CHARMM force-field in GROMACS format.
   - **Outputs**:
     1. `drug.itp`: Contains GROMACS itp.
     2. `drug.prm`: Parameters obtained from drug.str converted to GROMACS format and units.
     3. `drug.top`: A GROMACS topology file that incorporates (1) and (2).
     4. `drug_ini.pdb`: Coordinates of the molecule obtained from drug.mol2.
   - **Requirements**:
     - Python 3.5.2 or higher.
     - Numpy
     - Networkx (version 1.x series)
   - **Note**: Ensure to use the same version of CGenFF in your simulations that was used during parameter generation to avoid parameter mismatches.

2. **sort_mol2_bonds.pl**
   - **Description**: This Perl script reorders the bonds listed in a .mol2 file's @<TRIPOS>BOND section to follow specific conventions:
     1. Atoms on each line are in increasing order.
     2. Bonds appear in order of ascending atom number.
     3. Bonds involving the same atom in the first position appear in order of ascending second atom number.
   - **Usage**:
     ```sh
     perl sort_mol2_bonds.pl input.mol2 output.mol2
     ```
   - **Inputs**:
     1. `input.mol2`: The input .mol2 file.
     2. `output.mol2`: The output .mol2 file with sorted bonds.
   - **Outputs**: Sorted .mol2 file with bonds reordered according to the conventions mentioned above.
