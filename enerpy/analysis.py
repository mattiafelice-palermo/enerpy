import numpy as np
from functools import partial
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from .utils import convert_energy, normalize_data


def process_data_and_plot(parameter_id, parameter_value, cwd):
    """Processes and plots energy data from two files."""
    parameter_id = parameter_id.replace(' ', '_')
    energy_path = f"{cwd}/{parameter_id}/energy.xvg"
    qm_path = f"{cwd}/{parameter_id}/{parameter_id}.relaxscanact.dat"

    # Read and convert data
    mm_data = read_energy_data(energy_path)
    qm_data = read_energy_data(qm_path)
    mm_data[:, 1] = convert_energy(mm_data[:, 1], 'kJ/mol', 'kcal/mol')
    qm_data[:, 1] = convert_energy(qm_data[:, 1], 'Hartree', 'kcal/mol')

    # Normalize data
    mm_data = normalize_data(mm_data)
    qm_data = normalize_data(qm_data)

    def harmonic_spring(x, k, x0):
        """Harmonic spring model function."""
        return 0.5 * k * (x - x0)**2
    
    eq_value = float(parameter_value[0])
    
    if parameter_id.startswith('b'):
        eq_value *= 10 # Convert from nm to Angstrom
    
    harmonic_spring_fixed = partial(harmonic_spring, x0=eq_value)

    mm_k, mm_x0 = curve_fit(lambda x, k: harmonic_spring_fixed(x, k), qm_data[:, 0], mm_data[:, 1])
    qm_k, qm_x0 = curve_fit(lambda x, k: harmonic_spring_fixed(x, k), qm_data[:, 0], qm_data[:, 1])

    # Generate fitted data points for plotting
    x_fitted = np.linspace(min(qm_data[:, 0]), max(qm_data[:, 0]), 300)
    mm_fitted = harmonic_spring(x_fitted, *mm_k, eq_value)
    qm_fitted = harmonic_spring(x_fitted, *qm_k, eq_value)

    np.savetxt(f"{cwd}/{parameter_id}_mm_data.csv", np.column_stack((qm_data[:, 0], mm_data[:, 1])), delimiter=',', header='Bond Length,Energy', comments=f'{mm_k}, {mm_x0}')
    np.savetxt(f"{cwd}/{parameter_id}_qm_data.csv", qm_data, delimiter=',', header='Bond Length,Energy', comments=f'{qm_k}, {qm_x0}')

    # Plot data
    plt.figure(figsize=(10, 5))
    plt.plot(qm_data[:, 0], mm_data[:, 1], 'o-', label=f'MM (k={mm_k[0]:.2f} kcal/mol/Å²)')
    plt.plot(qm_data[:, 0], qm_data[:, 1], 's-', label=f'QM (k={qm_k[0]:.2f} kcal/mol/Å²)')
    plt.plot(x_fitted, mm_fitted, 'b--', label='Fitted MM')
    plt.plot(x_fitted, qm_fitted, 'r--', label='Fitted QM')
    plt.title(f'{parameter_id}')
    plt.xlabel('Bond Length (some unit)')
    plt.ylabel('Energy (kcal/mol)')
    plt.legend()
    plt.grid(True)

    plt.savefig(f"{cwd}/{parameter_id}_energy_plot.png")

    plt.close()

    # Extracting energy values from the data arrays
    mm_energies = mm_data[:, 1]
    qm_energies = qm_data[:, 1]

    cumulative_error = (mm_energies - qm_energies).sum()

    # Calculate MAE
    mae = np.mean(np.abs(mm_energies - qm_energies))

    # Calculate RMSE
    rmse = np.sqrt(np.mean((mm_energies - qm_energies)**2))

    # Calculate R^2
    ss_res = np.sum((mm_energies - qm_energies)**2)
    ss_tot = np.sum((qm_energies - np.mean(qm_energies))**2)
    r_squared = 1 - (ss_res / ss_tot)

    return {'mm_k': mm_k[0], 'qm_k': qm_k[0], 'mae': mae, 'rmse': rmse, 'r_squared': r_squared, 'cumulative': cumulative_error} 


def read_energy_data(file_path):
    """Reads energy data from a file, skipping headers."""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith(('#', '@')):
                parts = line.split()
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)
