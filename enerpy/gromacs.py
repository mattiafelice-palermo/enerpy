import os
import jobdispatcher as jd
import subprocess
import textwrap
import itertools

from iotools import mdp_to_dictionary, get_atom_coordinates
from utils import calculate_distance


def run_gromacs(parameters, topology, threads, total_cores):
    cwd = os.getcwd()

    def run_job(parameter_id, cwd, topology, threads):
        path = os.path.join(cwd, parameter_id)
        grompp(parameter_id, f'{cwd}/{topology}', path)
        mdrun(parameter_id, path, threads)
        extract_energy(parameter_id, path)      


    write_gromacs_mdp(parameters, cwd)

    jobs = []

    for parameter in parameters:
        parameter_id = parameter.replace(' ', '_')
        job = jd.Job(name=f"{parameter_id}", function=run_job, arguments = [parameter_id, cwd, topology, threads], cores=threads)
        jobs.append(job)

    dispatcher = jd.JobDispatcher(jobs, maxcores=total_cores, engine="multiprocessing")
    results = dispatcher.run()

def write_gromacs_mdp(parameters, cwd):
    mdp_file_content = textwrap.dedent(
        f""";--- RUN CONTROL
            integrator              = md-vv
            nsteps                  = 10000000 ; IGNORED
            dt                      = 0.0001
            ;--- OUTPUT CONTROL
            nstvout                 = 10000
            nstenergy               = 10000
            nstcalcenergy           = 1
            nstlog                  = 10000
            nstcomm                 = 1
            ;NEIGHBOR SEARCHING
            cutoff-scheme           = verlet
            pbc                     = xyz
            ;--- ELECTROSTATICS
            coulombtype             = PME
            rcoulomb                = 1.1 ; uncomment if you're not going to use mdrun -tunepme
            ;--- VAN DER WAALS
            vdwtype                 = cut-off
            rvdw                    = 1.1
            vdw-modifier            = potential-switch
            rvdw-switch             = 1.05
            ;--- BONDS
            constraints             = h-bonds
            constraint-algorithm    = lincs
        """)

    with open(f'{cwd}/resample.mdp', 'w') as f:
        f.write(mdp_file_content)

    for parameter in parameters:
        parameter_id = parameter.replace(' ', '_')  # Convert parameter into a valid directory/filename
        link_path = os.path.join(cwd, parameter_id, 'resample.mdp')  # Define full path for the new symbolic link

        # Check if the link already exists, if so, remove it first
        if os.path.islink(link_path):
            os.unlink(link_path)
        
        # Create a symbolic link pointing to 'target' with the name 'link_name' in the directory 'parameter_id'
        os.symlink(f'{cwd}/resample.mdp', link_path)

def grompp(parameter_id, topology, path, scan_step=None):
    if scan_step is None:
        scan_step = ""
    else:
        scan_step = f'_{scan_step}'
        topology_base =  os.path.basename(topology)[:-4]
        topology = f'{topology_base}{scan_step}.top'
    # Run GROMPP
    try:
        grompp = subprocess.run(f'module load gromacs && gmx grompp -f resample.mdp'
                                f' -c {parameter_id}{scan_step}.g96 -p {topology} -o resample{scan_step}.tpr -maxwarn 2 > grompp{scan_step}.out 2> grompp{scan_step}.err',
                                shell=True, cwd=path, check=True, executable='/bin/bash')
    except:
        raise RuntimeError(f"Failed to run GROMACS grompp command in {parameter_id}{scan_step}")

def mdrun(parameter_id, path, threads, scan_step=None):
    if scan_step is None:
        scan_step = ""
        rerun_flag = f'-rerun {parameter_id}.g96'
    else:
        scan_step = f'_{scan_step}'
        rerun_flag = ""
    # Run MDRUN
    try:
        mdrun = subprocess.run(f'module load gromacs && gmx mdrun -nt {threads}'
                               f' -s resample{scan_step}.tpr -deffnm {parameter_id}{scan_step} {rerun_flag} > mdrun{scan_step}.out 2> mdrun{scan_step}.err', shell=True, cwd=path, check=True, executable='/bin/bash')

    except:
        raise RuntimeError(f"Failed to run GROMACS mdrun command in {parameter_id}{scan_step}")

def extract_energy(parameter_id, path, scan_step=None):
    if scan_step is None:
        scan_step = ""
    else:
        scan_step = f'_{scan_step}'
    try:
        result = subprocess.run(f'module load gromacs && gmx energy -f {parameter_id}{scan_step}.edr -o {parameter_id}{scan_step}.xvg', 
                                input='Potential\n\n', check=True, text=True, capture_output=True, 
                                cwd=path, executable='/bin/bash', shell=True)
    except Exception as e:
        print(e)
        raise RuntimeError(f"Failed to run GROMACS energy command in {parameter_id}{scan_step}")
    
    print(f"Energy extracted from {parameter_id}{scan_step}.edr")

def run_gromacs_relaxed(parameters, topology, threads, total_cores, scans, include_hydrogens, hydrogen_sn):
    cwd = os.getcwd()

    def run_job(parameter_id, cwd, topology, threads, scan_step=None):
        path = os.path.join(cwd, parameter_id)
        grompp(parameter_id, f'{cwd}/{topology}', path, scan_step)
        mdrun(parameter_id, path, threads, scan_step)
        extract_energy(parameter_id, path, scan_step)    


    write_gromacs_mdp_minimization(cwd, parameters)
    all_parameters, hydrogen_sn = mdp_to_dictionary(topology, include_hydrogens=include_hydrogens)

    jobs = []

    for parameter in parameters:
        for scan_step in range(scans):
            parameter_id = parameter.replace(' ', '_') 
            write_constrained_topology(parameter, scan_step, topology, cwd, include_hydrogens, hydrogen_sn, all_parameters)
            job = jd.Job(name=f"{parameter_id}_{scan_step}", function=run_job, arguments = [parameter_id, cwd, topology, threads], keyword_arguments = {'scan_step': scan_step}, cores=threads)
            jobs.append(job)

    dispatcher = jd.JobDispatcher(jobs, maxcores=total_cores, engine="multiprocessing")
    results = dispatcher.run()
    for parameter in parameters:
        append_energies(parameter, scans, cwd)

def append_energies(parameter, scans, cwd):
    path = os.path.join(cwd, parameter)
    parameter_id = parameter.replace(' ', '_') 

    xvg_lines = []

    for scan_step in range(scans):
        with open(f'{cwd}/{parameter_id}/{parameter_id}_{scan_step}.xvg', 'rb') as f:
            f.seek(-2, 2)  # Jump to the second last byte.
            while f.read(1) != b'\n':  # Keep reading backwards until you find the new line.
                f.seek(-2, 1)
            xvg_lines.append(f.readline().decode()) # the last line of the file

    with open(f'{cwd}/{parameter_id}/energy.xvg', 'w') as f:
        for line in xvg_lines:
            f.write(line)

def write_constrained_topology(parameter_id, scan_step, topology, cwd, include_hydrogens, hydrogen_sn, parameters):
    # 1 calculate bond(s) distances
    if parameter_id.startswith('b'):
        atoms = [int(parameter_id.split()[1]), int(parameter_id.split()[2])]
    else:
        atoms = [int(parameter_id.split()[1]), int(parameter_id.split()[2]), int(parameter_id.split()[3])] 

    bonded_atoms = [ [int(key.split()[1]), int(key.split()[2])] for key in list(parameters.keys()) if key.startswith('b') ]

    # Check if any of the atoms in the bonded interaction are hydrogens 
    constrained_hydrogens = []
    for atom in atoms:
        if atom in hydrogen_sn:
            constrained_hydrogens.append(atom)

    if len(constrained_hydrogens) > 0 and not include_hydrogens:
        atoms = []

    # if len(constrained_hydrogens) == 1 and include_hydrogens and parameter_id.startswith('a'):
    #     # find the other bonded atom
    #     hydrogen = constrained_hydrogens[0]
    #     atoms.remove(hydrogen)
    #     for bond in bonded_atoms:
    #         print(bond)
    #         if hydrogen in bond and atoms[0] in bond:
    #             atoms.remove(atoms[0])
    #             break
    #         if hydrogen in bond and atoms[1] in bond:
    #             atoms.remove(atoms[0])
    #             break
    #     atoms = (hydrogen, atoms[0])
    
    # if len(constrained_hydrogens) == 2 and include_hydrogens and parameter_id.startswith('a'):
    #     atoms = (constrained_hydrogens[0], constrained_hydrogens[1])

            

    parameter_id = parameter_id.replace(' ', '_')

    coordinates = get_atom_coordinates(f'{cwd}/{parameter_id}/{parameter_id}_{scan_step}.g96', atoms)

    distance_dict = {}
    # Generate all combinations of atom pairs from the list
    for atom1, atom2 in itertools.combinations(atoms, 2):
        # Retrieve coordinates for both atoms
        coord1 = coordinates[atom1]
        coord2 = coordinates[atom2]
        # Calculate distance and store it in the dictionary
        distance_dict[(atom1, atom2)] = calculate_distance(coord1, coord2)

    topology_root = topology[:-4]
    
    # 2 copy the topology file from cwd to the parameter folder
    os.system(f'cp {topology} {cwd}/{parameter_id}/{topology_root}_{scan_step}.top')

    # 3 add them at the end of the file, it works!
    with open(f'{cwd}/{parameter_id}/{topology_root}_{scan_step}.top', 'a') as f:
        if atoms:
            f.write('\n\n[ constraints ]\n')
        for atom1, atom2 in distance_dict:
            f.write(f'{atom1} {atom2} 1 {distance_dict[(atom1, atom2)]}\n')

def write_gromacs_mdp_minimization(cwd, parameters):
    mdp_file_content = textwrap.dedent(
        f""";--- RUN CONTROL
            integrator               = steep    ; Algorithm (steep = steepest descent minimization)
            emtol                    = 1000.0   ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
            emstep                   = 0.01     ; Minimization step size
            nsteps                   = 50000    ; Maximum number of (minimization) steps to perform
            ;--- OUTPUT CONTROL
            ;nstvout                 = 10000
            ;nstenergy               = 10000
            ;nstcalcenergy           = 1
            ;nstlog                  = 10000
            ;nstcomm                 = 1
            ;NEIGHBOR SEARCHING
            cutoff-scheme           = verlet
            pbc                     = xyz
            ;--- ELECTROSTATICS
            coulombtype             = PME
            rcoulomb                = 1.1 ; uncomment if you're not going to use mdrun -tunepme
            ;--- VAN DER WAALS
            vdwtype                 = cut-off
            rvdw                    = 1.1
            vdw-modifier            = potential-switch
            rvdw-switch             = 1.05
            ;--- BONDS
            ;constraints             = h-bonds
            ;constraint-algorithm    = lincs
        """)
    
    with open(f'{cwd}/minimize.mdp', 'w') as f:
        f.write(mdp_file_content)

    for parameter in parameters:
        parameter_id = parameter.replace(' ', '_')  # Convert parameter into a valid directory/filename
        link_path = os.path.join(cwd, parameter_id, 'resample.mdp')  # Define full path for the new symbolic link

        # Check if the link already exists, if so, remove it first
        if os.path.islink(link_path):
            os.unlink(link_path)
        
        # Create a symbolic link pointing to 'target' with the name 'link_name' in the directory 'parameter_id'
        os.symlink(f'{cwd}/minimize.mdp', link_path)
