import os
import subprocess
import jobdispatcher as jd

def write_orca_input(parameters, geometry, method=None, basis_set=None, charge=None, multiplicity=None, threads_per_calc= None, cpus=None, memory=None, scan_steps = None):
    """
    Writes an ORCA input file for each parameter set.
    """
    if method is None:
        method = 'b97-3c'
    if basis_set is None:
        basis_set = ''
    if charge is None:
        charge = 0
    if multiplicity is None:
        multiplicity = 1
    if threads_per_calc is None:    
        threads_per_calc = 1
    if cpus is None:
        cpus = 1
    if memory is None:
        memory = 4000

    cwd = os.getcwd()

    for parameter in parameters:
        # Create a folder for each parameter
        folder_name = parameter.replace(' ', '_')
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        splitted = [str(int(atom_id)-1) for atom_id in parameter.split()[1:]]
        scan_type = parameter.split()[0].upper()
        atom_ids = ' '.join(splitted)
        eq_value = float(parameters[parameter][0][0])

        if scan_type == 'B':
            eq_value *= 10
        start = eq_value - eq_value * 0.1
        end = eq_value + eq_value * 0.1
        steps = scan_steps

        input_path = os.path.join(cwd, folder_name, f'{folder_name}.inp')
        geometry_path = os.path.join(cwd, geometry)
        orca_input =  f"! {method} {basis_set} Opt\n%pal nprocs {threads_per_calc} end\n%maxcore {memory}\n"
        orca_input += f"%geom Scan\n{scan_type} {atom_ids} = {start}, {end}, {steps}\nend\nend\n"
        orca_input += f"* xyzfile {charge} {multiplicity} {geometry_path}\n"
        with open(input_path, 'w') as f:
            f.write(orca_input)

def launch_orca_jobs(parameters, threads, total_cores, modulename=None):

    cwd = os.getcwd()
    def run_orca(parameter_id):
        folder_path = os.path.join(cwd, parameter_id)
        # Launch the job
        if modulename is None:
            module_command = ""
        else:
            module_command = f"module load {modulename} 2>/dev/null || true; "
        try:
            # Run the ORCA calculation, skip module load if not needed
            orca_command = f"{module_command} $(which orca) {parameter_id}.inp > {parameter_id}.out 2> {parameter_id}.err"
            subprocess.run(orca_command, shell=True, check=True, cwd=folder_path, executable='/bin/bash')
            print("ORCA job completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to complete ORCA job: {e}")
        
    jobs = []
    for parameter in parameters:
        parameter_id = parameter.replace(' ', '_')
        job = jd.Job(name=f"{parameter_id}", function=lambda parameter_id=parameter_id : run_orca(f"{parameter_id}"), cores=threads)
        jobs.append(job)

    dispatcher = jd.JobDispatcher(jobs, maxcores=total_cores, engine="multiprocessing")
    results = dispatcher.run()

    return results