import re
amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    takes a single peptide sequence, summs up the amino-acid masses and returns a dict of peptide and its according mass

    Parameters
    ----------
    peptide_seq: str
        contains peptide sequence
    amino_acid_mass_dict: dict
        contains amino-acid masses of all amino-acids

    Returns
    -------
    dict of {string, float}
        {pep, mass}

    """
    mass = sum(amino_acid_mass_dict[a] for a in peptide_seq)
    
    return {peptide_seq: mass}

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    takes a peptide sequence collection, summs up the amino-acid masses and returns a dict of peptides and their according masses

    Parameters
    ----------
    peptides: list
        contains peptide sequences
    amino_acid_mass_dict: dict
        contains amino-acid masses of all amino-acids

    Returns
    -------
    dict of {string, float}
        {pep, mass}
    """
    peptide_mass_collection = {}

    for pep in peptides:
        mass = 0

        for a in pep:
            mass = sum(amino_acid_mass_dict[a] for a in pep)
    
            peptide_mass_collection[pep] = mass

    return peptide_mass_collection

def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    mapps peptides to their m/z values

    Parameters
    ----------
    peptide_mass_map: dict
        contains peptides ansd their according mass
    charge : int
        Default charge state (default=2)
    proton_mass : float
        Mass of a proton (default=1.007 Da)
    Returns
    -------
    dict of {string, float}
        {pep, m/z}
    """
    mz_map = {}

    for pep, mass in peptide_mass_map.items():
        mz = (mass + charge * proton_mass) / charge
        mz_map[pep] = mz

    return mz_map


import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42, title = "ms spectrum"):
    """
    plots bar chart of m/z values and optional random_count_range

    Parameters
    ----------
    mz_values: list of float
        list of m/z values of peptides
    random_count_range: tuple of (int, int), optional
        Min and max of random intensities for the bars (default=(0,30000)).
    seed : int, optional
        Random seed for reproducibility (default=42).
    Returns
    -------
    bar chart
    """
    np.random.seed(seed)

    n_peaks = len(mz_values)

    intensities = np.random.randint(random_count_range[0],
                                    random_count_range[1],
                                    size=n_peaks)

    plt.figure(figsize=(10,5))
    plt.bar(mz_values, intensities, width=1, color='lavender')

    plt.xlabel('m/z')
    plt.ylabel('random Intensity')
    plt.title(title)
    plt.grid(axis='y', alpha=0.7)

    plt.show()

def fragment_peptide(peptide):
    """
    returns list of b- and y-ion sequences for any peptide input

    Parameters
    ----------
    peptide: string
        peptide sequence
    Returns
    -------
    list:
        b- and y-ion sequences
        
    """

    b_ions = [peptide[:i] for i in range(1, len(peptide))]

    y_ions = [peptide[-i:] for i in range(1, len(peptide))]

    fragment_list_peptide = b_ions + y_ions

    return(fragment_list_peptide)
