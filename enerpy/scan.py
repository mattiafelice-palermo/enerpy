import os
import argparse
from tabulate import tabulate

from iotools import mdp_to_dictionary, xyz_to_g96
from orca import write_orca_input, launch_orca_jobs
from gromacs import run_gromacs, run_gromacs_relaxed

from analysis import process_data_and_plot

def scan_workflow(args):
    # Read parameters and store in a dictionary
    parameters, hydrogen_sn = mdp_to_dictionary(args.topology, include_hydrogens=args.include_hydrogens)

    if not args.skip_orca_calc:
        # Create a ORCA file for each parameter
        write_orca_input(parameters, args.input_file, threads_per_calc=args.threads_per_calc, cpus=args.total_cores, scan_steps=args.scan_steps, charge=args.charge, multiplicity=args.multiplicity)   

        # Create jobs, each cding in the directory and calling orca input and launch through jobdispatcher
        launch_orca_jobs(parameters, args.threads_per_calc, args.total_cores)

    # Convert xyz files to g96 to keep fine geometry details
    

    if args.constrained_opt:
        # Convert xyz files to g96 to keep fine geometry details - all geometries in a single g96 file
        # TODO disable lincs constraint by default
        xyz_to_g96(parameters)
        # Option 1: recalculate the energies with MM force field
        run_gromacs(parameters, args.topology, args.threads_per_calc, args.total_cores)
    # Option 2: for each geometry, constrain bonded interaction and minimize
    else:
        xyz_to_g96(parameters, split=True)
        run_gromacs_relaxed(parameters, args.topology, args.threads_per_calc, args.total_cores, args.scan_steps, args.include_hydrogens, hydrogen_sn)

    # Minimize difference to obtain optimized parameters
    
    # Process data
    results = {}
    for parameter in parameters:
        result = process_data_and_plot(parameter, parameters[parameter][0], os.getcwd())
        results[parameter] = result

    # Prepare data for tabulation
    headers = ['Parameter', 'MM K (kcal/mol/Å²)', 'QM K (kcal/mol/Å²)', 'Cumulative Error (kcal/mol)', 'MAE (kcal/mol)', 'RMSE (kcal/mol)', 'R²']
    table = []

    for parameter, metrics in results.items():
        # Create a row for each parameter, extracting each metric
        row = [parameter, metrics['mm_k'], metrics['qm_k'], metrics['cumulative'], metrics['mae'], metrics['rmse'], metrics['r_squared']]
        table.append(row)

    # Sort the table by the MAE column (index 3)
    table = sorted(table, key=lambda x: x[3])

    # Print the table
    print(tabulate(table, headers=headers))#, tablefmt='grid'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scans a set of parameters and creates a folder for each parameter.')
    parser.add_argument('-i', '--input-file', type=str, help='Structure from which to start the calculations from.')
    parser.add_argument('-t', '--topology', type=str, help='Path to the input file containing the parameters to scan.')
    parser.add_argument('-tpc', '--threads-per-calc', type=int, help='Number of threads to use for each calculation.')
    parser.add_argument('--total-cores', type=int, help='Number of cores to use for the calculations.')
    parser.add_argument('--scan-steps', type=int, default=8, help='Number of steps to scan (default: 8).')
    parser.add_argument('--skip-orca-calc', action='store_true', help='Skip the ORCA calculation.')
    parser.add_argument('--include-hydrogens', action='store_true', help='Include hydrogens in the calculation.')
    parser.add_argument('--constrained-opt', action='store_true', help='Perform a constrained optimization.')
    parser.add_argument('--charge', type=int, help='Charge of the simulated specie.')
    parser.add_argument('--multiplicity', type=int, help='Multiplicity of the simulate specie.')
    args = parser.parse_args()
    scan_workflow(args)
