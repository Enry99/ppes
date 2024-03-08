#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
@ Program: ppes
@ Author: Enrico Pedretti
@ Created: 2024-04-03
'''

from ase.io import read
import numpy as np
import os, sys, shutil, json
from ase.calculators.espresso import Espresso
from ase.calculators.vasp import Vasp
from ase.io.espresso import read_fortran_namelist
from ase.io.vasp import write_vasp
import ase_custom

TEST = False

# Define the folder names for the isolated structures
upper_folder = 'isolated/upper'
lower_folder = 'isolated/lower'
sweep_folder = 'distance_sweep'


def write_input(atoms, folder, settings_dict, filename_label=None):
    """
    Write input files for different electronic structure calculation programs.

    Parameters:
        atoms (ase.Atoms): The atomic structure.
        folder (str): The folder path where the input files will be written.
        settings_dict (dict): A dictionary containing the settings for the electronic structure calculation.
        filename_label (str, optional): The label for the input file. Defaults to None.

    Returns:
        None
    """

    folder = folder.rstrip('/') #remove possible trailing slash in folder path
    

    if settings_dict["program"] == 'qe':
        dft_settings_dict = parse_espresso_settings(settings_dict["pwi_head"])

        calculator = Espresso(label = filename_label if filename_label else 'input',
                            directory=folder,
                            pseudopotentials=dft_settings_dict['pseudopotentials'],
                            kpts=dft_settings_dict["kpts"],
                            koffset=dft_settings_dict["koffset"],
                            input_data=dft_settings_dict)
        calculator.write_input(atoms)
    

    elif settings_dict["program"] == 'vasp':
        
        if settings_dict.get("poscar_only", False):
            if not os.path.exists(folder):
                os.makedirs(folder, exist_ok=True)
            write_vasp(f'{folder}/POSCAR', atoms, sort=False)
            shutil.copyfile('INCAR', f'{folder}/INCAR')
            shutil.copyfile('KPOINTS', f'{folder}/KPOINTS')
            shutil.copyfile('POTCAR', f'{folder}/POTCAR')

        else:
            os.environ["VASP_PP_PATH"] = settings_dict["vasp_pp_path"] #set the path to the VASP pseudopotentials
            
            calculator = Vasp(atoms=atoms, 
                            directory=folder,
                            xc='PBE', 
                            setups=settings_dict["vasp_pseudo_setups"])
            calculator.read_incar('INCAR')
            if 'magmom' in calculator.list_float_params:
                if len(calculator.list_float_params['magmom']) != len(atoms):
                    calculator.list_float_params['magmom'] = [0.5]*len(atoms) #TODO: make this a setting
            calculator.read_kpoints('KPOINTS')

            calculator.write_input(atoms)


def parse_espresso_settings(filename: str):
    """
    Parses the espresso settings from a given file and returns a dictionary containing the parsed settings.

    Parameters:
    filename (str): The path to the file containing the espresso settings.

    Returns:
    dict: A dictionary containing the parsed espresso settings.

    """

    with open(filename, 'r') as f:
        block_str_list = f.readlines()


    # parse namelist section and extract remaining lines
    dftprogram_settings_dict, card_lines = read_fortran_namelist(block_str_list)

    # parse ATOMIC_SPECIES and K_POINTS
    for i, line in enumerate(card_lines):
        if 'ATOMIC_SPECIES' in line.upper():
            atomic_species_index = i
        if 'K_POINTS' in line.upper():
            k_points_index = i

    # ATOMIC_SPECIES
    i = atomic_species_index + 1

    dftprogram_settings_dict['pseudopotentials'] = {}
    while i < len(card_lines):
        line = card_lines[i]

        if ('ATOMIC_POSITIONS' in line.upper() or
                'K_POINTS' in line.upper() or
                'ADDITIONAL_K_POINTS' in line.upper() or
                'CELL_PARAMETERS' in line.upper() or
                'CONSTRAINTS' in line.upper() or
                'OCCUPATIONS' in line.upper() or
                'ATOMIC_VELOCITIES' in line.upper() or
                'ATOMIC_FORCES' in line.upper() or
                'SOLVENTS' in line.upper() or
                'HUBBARD' in line.upper()):
            break

        element, mass, pseudo = line.split()
        dftprogram_settings_dict['pseudopotentials'].update({element: pseudo})
        i += 1

    # K_POINTS
    if 'gamma' in card_lines[k_points_index].split()[1].strip().lower():
        dftprogram_settings_dict['kpts'] = None
        dftprogram_settings_dict['koffset'] = None
    else:
        line = card_lines[k_points_index + 1]
        dftprogram_settings_dict['kpts'] = list(map(int, line.split()[:3]))
        dftprogram_settings_dict['koffset'] = list(map(int, line.split()[3:]))

    return dftprogram_settings_dict


def launch_job(folder, program, jobscript_file, input_label=None):
    """
    Launches a job by copying a jobscript file to a specified folder and executing it.

    Args:
        folder (str): The path to the folder where the jobscript file will be copied.
        program (str): The program to be used for executing the jobscript. Supported values are 'qe' and 'vasp'.
        jobscript_file (str): The path to the jobscript file.
        input_label (str, optional): An optional label to be used in the jobscript execution command. Defaults to None.

    Returns:
        None
    """

    folder = folder.rstrip('/') #remove possible trailing slash in folder path
    
    shutil.copyfile(jobscript_file, f'{folder}/jobscript')
    main_folder = os.getcwd()

    os.chdir(folder)
    if program == 'qe':
        if TEST:
            print(f'sbatch jobscript {input_label}.pwi {input_label}.pwo from {folder}')
        else:
            os.system(f'sbatch jobscript {input_label}.pwi {input_label}.pwo')
    elif program == 'vasp':
        if TEST:
            print(f'sbatch jobscript from {folder}')
        else:
            os.system(f'sbatch jobscript')
    os.chdir(main_folder)


def sweep_coordinate(settings_dict):
    """
    Determines the coordinate to sweep based on the given settings dictionary.

    Args:
        settings_dict (dict): A dictionary containing the settings.

    Returns:
        int: The index of the coordinate to sweep.

    """
    which_coord = settings_dict.get("which_coord", "z")
    if which_coord == 'x':
        which_coord = 0
    elif which_coord == 'y':
        which_coord = 1
    elif which_coord == 'z':
        which_coord = 2

    return which_coord


def distance(settings_dict, upper, lower):
    """
    Calculate the distance between the upper and lower atoms.

    Args:
        upper (Atoms): The upper atoms object.
        lower (Atoms): The lower atoms object.

    Returns:
        float: The distance between the upper and lower atoms.
    """

    which_coord = sweep_coordinate(settings_dict)

    
    if "upper_range_fordistance" in settings_dict or "upper_species_fordistance" in settings_dict or "upper_thickness_fordistance" in settings_dict:
        #These options are combined (cascade) in this order

        if "upper_range_fordistance" in settings_dict:
            upper_indices = [i for i in range(len(upper)) if upper.positions[i,which_coord] > settings_dict["upper_range_fordistance"][0]  \
                            and upper.positions[i,which_coord] < settings_dict["upper_range_fordistance"][1]]
            upper = upper[upper_indices]
        
        if "upper_species_fordistance" in settings_dict:
            upper_indices = [atom.index for atom in upper if (atom.symbol in settings_dict["upper_species_fordistance"])]
            upper = upper[upper_indices]

        if "upper_thickness_fordistance" in settings_dict:
            zmin = np.min(upper.positions[:,which_coord])
            upper_indices = [i for i in range(len(upper)) if upper.positions[i,which_coord] < zmin + settings_dict["upper_thickness_fordistance"]]
            upper = upper[upper_indices]

        z_upper = np.mean(upper.positions[:,which_coord])

    else: z_upper =  np.min(upper.positions[:,which_coord])


    if "lower_range_fordistance" in settings_dict or "lower_species_fordistance" in settings_dict:
        #These options are combined (cascade) in this order

        if "lower_range_fordistance" in settings_dict:
            lower_indices = [i for i in range(len(lower)) if lower.positions[i,which_coord] > settings_dict["lower_range_fordistance"][0]  \
                            and lower.positions[i,which_coord] < settings_dict["lower_range_fordistance"][1]]
            lower = lower[lower_indices]
        
        if "lower_species_fordistance" in settings_dict:
            lower_indices = [atom.index for atom in lower if (atom.symbol in settings_dict["lower_species_fordistance"])]
            lower = lower[lower_indices]

        if "lower_thickness_fordistance" in settings_dict:
            zmax = np.max(lower.positions[:,which_coord])
            lower_indices = [i for i in range(len(lower)) if lower.positions[i,which_coord] > zmax - settings_dict["lower_thickness_fordistance"]]
            lower = lower[lower_indices]

        z_lower = np.mean(lower.positions[:,which_coord])
    
    else: z_lower =  np.max(lower.positions[:,which_coord])

    return z_upper - z_lower


def generate_isolated(settings_dict):
    """
    Generate isolated upper and lower structures based on the provided settings.

    Args:
        settings_dict (dict): A dictionary containing the settings for generating the isolated structures.
    """

    if "filename_upper" in settings_dict and "filename_lower" in settings_dict:
        # Upper and lower already present, no need for selection
        
        upper = read(settings_dict["filename_upper"])
        lower = read(settings_dict["filename_lower"])

        # If both upper and lower already calculated, just return
        if hasattr(upper, "calc") and hasattr(lower, "calc") and upper.calc is not None and lower.calc is not None:
            if upper.calc.results.get("energy", None) is not None and lower.calc.results.get("energy", None) is not None:
                return upper, lower, None, None
            
    else:
        # Upper and lower not present, need to select from full structure
        atoms = read(settings_dict["filename_full"])

        if "upper_atoms" in settings_dict:
            upper_indices = settings_dict["upper_atoms"]
        elif "upper_species" in settings_dict:
            upper_indices = [atom.index for atom in atoms if (atom.symbol in settings_dict["upper_species"])]
        elif "upper_min" in settings_dict:
            which_coord = sweep_coordinate(settings_dict)
            upper_indices = [i for i in range(len(atoms)) if atoms.positions[i,which_coord] > settings_dict["upper_min"]]
        else:
            raise RuntimeError("No method to select upper atoms provided.")

        lower_indices = [i for i in range(len(atoms)) if i not in upper_indices]

        upper = atoms[upper_indices]
        lower = atoms[lower_indices]
    

    write_input(upper, upper_folder, settings_dict, "upper")
    write_input(lower, lower_folder, settings_dict, "lower")

    if settings_dict.get("launch_isolated", False):
        if settings_dict.get("poscar_only", False):
            raise ValueError("Cannot launch isolated structures with poscar_only set to True.")
        launch_job(upper_folder, settings_dict["program"], settings_dict["jobscript_file"], "upper")
        launch_job(lower_folder, settings_dict["program"], settings_dict["jobscript_file"], "lower")
    else:
        print("Upper and lower were written to isolated/upper and isolated/lower, but NOT calculated, as requested.")

    if settings_dict.get("keep_order", settings_dict.get("poscar_only", False)):
        return upper, lower, upper_indices, atoms
    else:
        return upper, lower, None, None
            

def generate_sweep(upper, lower, upper_indices, atoms, settings_dict):
    """
    Generate a sweep of atomic positions by moving the upper atoms along a specified coordinate.

    Args:
        upper (Atoms): The upper atoms object.
        lower (Atoms): The lower atoms object.
        settings_dict (dict): A dictionary containing the settings for the sweep.

    Returns:
        None
    """

    range_sweep = settings_dict["range_sweep"]
    step = settings_dict["step"]

    which_coord = sweep_coordinate(settings_dict)

    initial_distance = distance(settings_dict, upper, lower)

    if settings_dict.get("reference_distance", None) is not None:
        # If the initial distance is set, we need to move the upper atoms to the correct position
        transl = settings_dict["reference_distance"] - initial_distance
        if atoms is None:  
            upper.positions[:,which_coord] += transl
        else:
            atoms.positions[upper_indices,which_coord] += transl
        
        initial_distance = settings_dict["reference_distance"]


    # sweep the upper atoms, keeping the lower atoms fixed
        
    z_values = []
    for dz in np.arange(range_sweep[0], range_sweep[1], step):
        if atoms is None:
            new_upper = upper.copy()
            new_upper.positions[:,which_coord] += dz
            new_atoms = new_upper + lower
        else:
            new_atoms = atoms.copy()
            new_atoms.positions[upper_indices,which_coord] += dz
        
        folder = "{0}/d_{1:4.2f}".format(sweep_folder, dz)
        os.makedirs(folder, exist_ok=True)
        write_input(new_atoms, folder, settings_dict, "d_{:4.2f}".format(dz))
        z_values.append(dz+initial_distance)
        
        if settings_dict.get("launch_sweep", False):
            launch_job(folder, settings_dict["program"], settings_dict["jobscript_file"], "d_{:4.2f}".format(dz))


    with open('distance.txt', 'w') as f:
        f.write('\n'.join(map(str, z_values)))


def get_results(settings_dict):
    """
    Get the results from the isolated and sweep calculations.

    Args:
        settings_dict (dict): A dictionary containing the settings for the electronic structure calculation.

    Returns:
        dict: A dictionary containing the results of the electronic structure calculation.
    """

    results = {}


    try:
        if "filename_upper" in settings_dict and "filename_lower" in settings_dict:  
            upper = read(settings_dict["filename_upper"])
            lower = read(settings_dict["filename_lower"])

            if hasattr(upper, "calc") and hasattr(lower, "calc"):
                if upper.calc.results.get("energy", None) is not None and lower.calc.results.get("energy", None) is not None:
                    results["upper_energy"] = upper.get_potential_energy()
                    results["lower_energy"] = lower.get_potential_energy()

        if "upper_energy" not in results or "lower_energy" not in results:
            up_fn = 'OUTCAR' if settings_dict["program"] == 'vasp' else 'upper.pwo'
            dw_fn = 'OUTCAR' if settings_dict["program"] == 'vasp' else 'lower.pwo'
            upper = read(f'{upper_folder}/{up_fn}')
            lower = read(f'{lower_folder}/{dw_fn}')
            results["upper_energy"] = upper.get_potential_energy()
            results["lower_energy"] = lower.get_potential_energy()
    except:
        print("Could not read isolated energies. Total energies will be displayed.")
        results["upper_energy"] = 0
        results["lower_energy"] = 0


    range_sweep = settings_dict["range_sweep"]
    step = settings_dict["step"]

    #read z values from file
    with open('distance.txt', 'r') as f:
        z_values = list(map(float, f.readlines()))

    z_values_present = []
    energies = []
    for dz, z in zip(np.arange(range_sweep[0], range_sweep[1], step), z_values):
        folder = "{0}/d_{1:4.2f}".format(sweep_folder, dz)
        if settings_dict["program"] == 'qe':
            file = "{0}/d_{1:4.2f}.pwo".format(folder, dz)
        elif settings_dict["program"] == 'vasp':
            file = f'{folder}/OUTCAR'
    
        if os.path.exists(file):
            try:
                atoms = read(file)
                energy = atoms.get_potential_energy()
            except:
                continue

            z_values_present.append(z)
            energies.append(energy - results["upper_energy"] - results["lower_energy"])


    results["z_values"] = z_values_present
    results["energies"] = energies

    return results


def plot_results(results):
    """
    Plot the PPES

    Args:
        results (dict): A dictionary containing the PPES energies.
    """

    import matplotlib.pyplot as plt

    plt.plot(results["z_values"], results["energies"], '-o')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Energy (eV)')
    plt.savefig('ppes.png', dpi=300, bbox_inches='tight')




def main():
    
    with open('settings.json', 'r') as f:
        settings_dict = json.load(f)


    command = sys.argv[1]

    if command == 'gen':
        print("Generating isolated structures...")
        upper, lower, upper_indices, atoms = generate_isolated(settings_dict)
        print("Isolated structures generated.")

        print("Generating PPES files...")
        generate_sweep(upper, lower, upper_indices, atoms, settings_dict)
        print("PPES files generated.")

    elif command == 'plot':
        results = get_results(settings_dict)
        plot_results(results)
        print("PPES plot generated.")

    else:
        print("Command not recognized.")