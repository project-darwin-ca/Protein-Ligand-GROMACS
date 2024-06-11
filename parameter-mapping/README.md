
# Parameter Mapping

This folder contains Molecular Dynamics Parameter (MDP) files for setting up and running different stages of molecular dynamics simulations using GROMACS. Each MDP file specifies various simulation parameters and settings.

## Files

1. **EM.mdp**
   - **Description**: Energy Minimization
   - **Purpose**: To relax the system to remove any steric clashes or inappropriate geometry before running molecular dynamics simulations.
   - **Key Parameters**: 
     - Minimization method: Steepest Descent
     - Maximum number of steps: 50000
     - Convergence criterion: 1000.0 kJ/mol/nm

2. **NVT.mdp**
   - **Description**: Constant Number, Volume, and Temperature (NVT) Equilibration
   - **Purpose**: To equilibrate the system at a constant volume and temperature.
   - **Key Parameters**:
     - Temperature coupling: V-rescale thermostat
     - Reference temperature: 300 K
     - Duration: 100 ps

3. **NPT.mdp**
   - **Description**: Constant Number, Pressure, and Temperature (NPT) Equilibration
   - **Purpose**: To equilibrate the system at constant pressure and temperature.
   - **Key Parameters**:
     - Pressure coupling: Parrinello-Rahman barostat
     - Reference pressure: 1 bar
     - Duration: 100 ps

4. **MD.mdp**
   - **Description**: Molecular Dynamics Production Run
   - **Purpose**: To perform the production run for molecular dynamics simulation.
   - **Key Parameters**:
     - Duration: 1 ns or longer
     - Temperature coupling: V-rescale thermostat
     - Pressure coupling: Parrinello-Rahman barostat

5. **ions.mdp**
   - **Description**: Ion Addition
   - **Purpose**: To prepare the system for ion addition by setting the necessary parameters.
   - **Key Parameters**:
     - Number of ions to add: Specify based on system requirements

## Usage

Each of these MDP files can be used with the `gmx grompp` command in GROMACS to preprocess the simulation input files, followed by the `gmx mdrun` command to run the simulations.

### Example Workflow

1. **Energy Minimization**
   ```bash
   gmx grompp -f EM.mdp -c input.gro -p topol.top -o em.tpr
   gmx mdrun -v -deffnm em
   ```

2. **NVT Equilibration**
   ```bash
   gmx grompp -f NVT.mdp -c em.gro -p topol.top -o nvt.tpr
   gmx mdrun -v -deffnm nvt
   ```

3. **NPT Equilibration**
   ```bash
   gmx grompp -f NPT.mdp -c nvt.gro -p topol.top -o npt.tpr
   gmx mdrun -v -deffnm npt
   ```

4. **Production MD Run**
   ```bash
   gmx grompp -f MD.mdp -c npt.gro -p topol.top -o md.tpr
   gmx mdrun -v -deffnm md
   ```

## License
This project is licensed under the MIT License. See the [LICENSE](../LICENSE) file for details.

## Contact
For questions, bug reports, or help, please contact the repository maintainers.

