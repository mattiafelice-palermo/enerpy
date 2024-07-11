import os
import re

def xyz_to_g96(parameters, split=False, box_dimensions=[20.0, 20.0, 20.0]):
    for parameter in parameters:
        parameter_id = parameter.replace(' ', '_') 
        xyz_filename = f'{parameter_id}/{parameter_id}.allxyz'
        xyz_path = os.path.join(os.getcwd(), xyz_filename)
        g96_filename = f'{parameter_id}/{parameter_id}'
        positions = read_xyz_positions(xyz_path)

        # Remove existing file to avoid appending to old content, if present
        g96_filename = f'{g96_filename}.g96'
        split_filename = f'{g96_filename[:-4]}'
        if os.path.exists(f'{g96_filename}'):
            print(f'Removing existing file {g96_filename}')
            os.remove(f'{g96_filename}')
        for i, position in enumerate(positions):
            if split:
                g96_filename =f'{split_filename}_{i}.g96' # -4 to remove extension
                if os.path.exists(f'{g96_filename}'):
                    os.remove(f'{g96_filename}')
            
            write_g96(position[2:], g96_filename) # Start from the third line to skip headers
    return g96_filename

def read_xyz_positions(xyz_filename):
    with open(xyz_filename, 'r') as file:
        xyz_positions = [[]]
        counter = 0
        for line in file:
            # Prepare a new list for the next set of positions upon encountering '>'
            if line.startswith('>'):
                counter +=1
                xyz_positions.append([])
                continue
            xyz_positions[counter].append(line)
    return xyz_positions

def write_g96(positions, output_filename, box_size=(20.0, 20.0, 20.0)):
    with open(output_filename, 'a') as file:
        file.write("TITLE\n")
        file.write("Converted from XYZ format\n")
        file.write("END\n")
        file.write("POSITIONRED\n")
        for i, position in enumerate(positions):
            pos = position.split()[1:]
            # Format the coordinates with sufficient space and precision
            formatted_pos = "{:15.9f}{:15.9f}{:15.9f}".format(float(pos[0])/10., float(pos[1])/10., float(pos[2])/10.)
            file.write(formatted_pos + "\n")
        file.write("END\n")
        file.write("BOX\n")
        file.write("{:15.9f}{:15.9f}{:15.9f}\n".format(box_size[0], box_size[1], box_size[2]))
        file.write("END\n")

def get_atom_coordinates(filename, atom_numbers):
    """
    Read a file and return the coordinates of specified atoms.

    :param filename: Path to the file containing atom coordinates.
    :param atom_numbers: A list of atom numbers (1-based index) for which coordinates are requested.
    :return: A dictionary with atom numbers as keys and their coordinates as values (tuples).
    """
    coordinates = {}
    start_reading = False
    atom_index = 1  # Start counting from 1 since atom_numbers are 1-based index
    
    with open(filename, 'r') as file:
        for line in file:
            if 'POSITIONRED' in line:
                start_reading = True
                continue
            if 'END' in line and start_reading:
                break
            if start_reading:
                if atom_index in atom_numbers:
                    # Extract coordinates from the line, convert them to float, and store them
                    parts = line.strip().split()
                    coordinates[atom_index] = tuple(float(part) for part in parts)
                atom_index += 1

    return coordinates

def mdp_to_dictionary(file_path, include_hydrogens=False) -> dict:
    """
    Converts MDP file contents to a pandas DataFrame representing the simulation parameters.

    Args:
        file_path (str): Path to the MDP file.

    Returns:
        dict: Dictionary containing the parameters extracted from the MDP file.
    """

    bond_lines = []
    angle_lines = []
    read_bonds = False
    read_angles = False
    read_atomtypes = False
    read_atoms = False

    # to exclude bonds containing hydrogens from the scan
    hydrogen_at = [] # atomtypes corresponding to an atom number of 1
    hydrogen_sn = [] # serial number of hydrogen atoms in the structure

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(';') or not line.strip(): # skip empty lines or comments
                continue
            if "[ atomtypes ]" in line:
                read_atomtypes = True
                continue
            if read_atomtypes and line.startswith('['):
                read_atomtypes = False
            if read_atomtypes:
                atom_number = int(line.split()[1])
                if atom_number == 1:
                    atom_type = line.split()[0]
                    hydrogen_at.append(atom_type)
            if "[ atoms ]" in line:
                read_atoms = True
                continue
            if read_atoms and line.startswith('['):
                read_atoms = False
            if read_atoms:
                atom_type = line.split()[1]
                if atom_type in hydrogen_at:
                    hydrogen_sn.append(int(line.split()[0]))
            if "[ bonds ]" in line:
                read_bonds = True
                continue
            if read_bonds:
                bond_lines.append(line)
            if read_bonds and line.startswith('['):
                read_bonds = False
            if "[ angles ]" in line:
                read_angles = True
                continue
            if read_angles:
                angle_lines.append(line)
            if read_angles and line.startswith('['):
                read_angles = False

    # Regex pattern to match numbers and floating point numbers, ignoring optional comments
    pattern_bond = r'\s*(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)(\s*;.*)?'
    pattern_angle = r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)(\s*;.*)?'

    params = {}

    # Iterate over each line in the data list to extract bond parameters
    for line in bond_lines:
        match = re.match(pattern_bond, line)
        if match:
            # Extract all groups, but filter out the comment
            values = match.groups()[:-1]  # Ignore the last group which is the comment      
            if (not include_hydrogens) and (int(values[0]) in hydrogen_sn or int(values[1]) in hydrogen_sn):
                print(f"Skipping bond {values[0]}-{values[1]} as it contains a hydrogen atom")
                continue
            key = f"b {values[0]} {values[1]}"
            params.setdefault(key, []).append([values[3], values[4]])

    # Iterate over each line in the angle data list to extract bond parameters
    for line in angle_lines:
        match = re.match(pattern_angle, line)
        if match:
            # Extract all groups, but filter out the comment
            values = match.groups()[:-1]  # Ignore the last group which is the comment
            key = f"a {values[0]} {values[1]} {values[2]}"
            params.setdefault(key, []).append([values[4], values[5]])
    return params, hydrogen_sn
