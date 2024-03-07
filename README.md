# pPES

Script to generate and launch PPES calculations with VASP and Quantum Espresso

## Installation

### Download and setup

First, clone the repository into your local machine: 

`git clone https://github.com/Enry99/ppes.git`

Once the download is completed, go inside the downloaded folder and run: 

`bash install.sh`

to add the executable in the `PATH` variable. With this operation, you can launch the program from any folder. 

### Dependencies

The program requires the following Python libraries:
- [`ase`](https://wiki.fysik.dtu.dk/ase/index.html)
- [`numpy`](https://numpy.org/)
- [`matplotlib`](https://matplotlib.org/)

The installation of the dependencies can be done for example with `pip`, using the command:

`python3 -m pip install -r requirements.txt`



## Program use

To generate the PPES structures, use the command:  

`ppes gen`

To plot the PPES profile, use the command:

`ppes plot`

All the settings must be provided in a `settings.json` file (an example with all possible flags is provided inside this repository)

### Settings parameters


**PPES settings**:
- `range_sweep`: range of translation [dz_min, dz_max]. The upper structure will be moved between [z_initial + dz_min, z_initial+dz_max]. The values can be negative.
- `step`: step for the translations
- `which_coord`: coordinate along which the translation is performed (default 'z'). It can be also 'x' or 'y'. All the variables referred to 'upper' belong to the part of the structure which is moved, while 'lower' belong to the part that stays still.
- `reference_distance`: initial distance associated with no translation (dz=0). The distance is calculated according to the corresponding flags, if specified, otherwise it it simply the difference between min(upper) and max(lower)

**Input structure(s)**:\
You can either provide a single file containing the upper and lower part using `filename_full`, specifying the flags for the selection of the upper part, or provide two separate files (`filename_upper` and `filename_lower`) with upper and lower structures.

- `filename_full`: path of the file with the full structure (must be in a ASE-readable format)
- `filename_upper`: path of the file with the upper structure (must be in a ASE-readable format)
- `filename_lower`: path of the file with the lower structure (must be in a ASE-readable format)


- `upper_atoms`: selection of atom indices for the "upper" structure, i.e. the one that is moved
- `upper_species`: select the "upper" structure by species
- `upper_min`: select the "upper" structure by the minimum coordinate (e.g. all the atoms with z>zmin if which_coord is 'z')


**Distance calculation:**\
The following settings are used to select the atoms belonging to upper and to lower to calculate the distance. 
If specified, the distance will be calculated as
mean(z_[upper subset]) - mean(z_[lower subset]).
If not present, it will simply be min(z_upper) - max(z_lower)
These options can be used in conjunction, as the subsetting is done by applying in this order the conditions.
- `upper_range_fordistance`: range [z_min, z_max] to select the upper atoms to average in order to calculate the distance
- `upper_species_fordistance`: species of the upper slab to average in order to calculate the distance
- `upper_thickness_fordistance`: window [min(z_upper), min(z_upper) + upper_thickness_fordistance] to select the upper atoms in order to calculate the distance

- `lower_range_fordistance`: range [z_min, z_max] to select the lower atoms to average in order to calculate the distance
- `lower_species_fordistance`: species of the lower slab to average in order to calculate the distance
- `lower_thickness_fordistance`: window [max(z_lower) - lower_thickness_fordistance, max(z_lower)] to select the lower atoms in order to calculate the distance


**DFT program settings**:
- `program`: "qe" or "vasp"

For QE:
- `pwi_head`: path for initial part of a pwi file, containing the settings for the calculations (system, electrons, kpoints, etc...)


For VASP:
- `vasp_pp_path`: the directory must contain a folder named 'potpaw_PBE' with all the folders with the name of the elements, each containing its POTCAR
- `vasp_pseudo_setups`: version of the potentials, as a dictionary, e.g. {"base": "recommended", "Li": "_sv"}
- `poscar_only`: "write only the POSCAR

Note that 'vasp_pp_path' must contain folders with ASE nomenclature, e.g. 'potpaw_PBE' (containing all the folders with the POTCARs)


**Job submission settings:**
- `launch_sweep`: if false, generate only the input file without submitting the sweep (PPES) calculations.
- `launch_isolated`: if false, generate only the files without submitting the calculations for the two isolated structures (upper and lower). Note that if the two structures are provided as a file that contains the energy, e.g. OUTCAR. The calculation is not re-done and that value is used.
- `jobscript_file`. Path of a working jobscript file. For QE you need to include `-input $1 >> $2` in the line where pw.x is run.


