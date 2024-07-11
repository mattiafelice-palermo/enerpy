import math
import numpy as np

def convert_energy(values, from_unit, to_unit):
    """Converts energy from one unit to another."""
    if from_unit == 'kJ/mol' and to_unit == 'kcal/mol':
        return values * 0.239005736  # conversion factor
    elif from_unit == 'Hartree' and to_unit == 'kcal/mol':
        return values * 627.509  # conversion factor
    return values

def normalize_data(data):
    """Normalizes data by subtracting the minimum value."""
    min_value = np.min(data[:, 1])
    data[:, 1] -= min_value
    return data



def calculate_distance(p1, p2):
    """Calculate the Euclidean distance between two points in 3D."""
    return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2 + (p2[2] - p1[2]) ** 2)