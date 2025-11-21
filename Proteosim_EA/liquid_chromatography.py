import re
from pyteomics import achrom

def predict_lc_retention_times(peptides):
    """
    returns a dictionary which contains each peptide and its predicted retention time

    Parameters
    ----------
    peptides: list
        list of peptide sequences

    Returns
    -------
    dict of {str: float}
        dictionary mapping each peptide to its predicted retention time
        retention time rounded to two decimals
    """

    retentions_times_peptides = {}

    for pep in peptides:
        relative_RT = achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0)
        retentions_times_peptides[pep] = float(round(relative_RT, 2)) #auf zwei Nachkommastellen runden


    return retentions_times_peptides


import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Plots a histogram of retention times.

    Parameters
    
    retention_times: list
        List of retention times (floats or ints)
    resolution: int, optional
        Number of bins in the histogram (default is 30)

    Returns
    histogram

    """
    plt.figure() # das vorherhige Histogramme nicht überschrieben werden
    #Histogramm erstellen
    plt.hist(retention_times, bins=resolution, color='lavender', edgecolor='lavender')

    #Achsenbeschriftungen
    plt.xlabel("Retention Time")
    plt.ylabel("Frequency")

    #Titel
    plt.title("Retention Times")

    #Gitterlinien
    plt.grid(axis='y', alpha=0.75)

    plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    sets a retention-time window and returns only the peptides that fall inside that window

    Parameters
    peptide_rt_map: dict
        contains the retention times of the peptides and the peptides
    lower_ret_time: float
        sets the lower border
    upper_ret_time: float
        sets the upper border

    Returns
    peptides_in_window: list
        peptides that fall inside that window

    """
    peptides_in_window = {}
    
    for pep, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time:
            peptides_in_window[pep] = rt

    

    return peptides_in_window