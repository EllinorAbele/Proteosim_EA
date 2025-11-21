import re
enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)', # nach Lysin geschnitten
    'LysN': r'(?=K)',  # vor Lysin geschnitten
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',       # eckige Klammer ist entweder oder
# schneidet am C-Terminus von Lysin ODER Arginin, und Befehl noch weiter "vor" zu schauen, nur wenn da kein Prolin drauf folgt
}

def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Digest a protein sequence according to a specified enzyme & cleavage pattern
    Filter peptide length

    Parameters
    protein_seq : str
        Protein sequence to be digested
        ....
    cleave_pattern
    min_pep_len, default 5
    max_pep_len, default 30

    Returns
    filtered list of peptides

    """
    peptides =[]
    filtered_peptides = []

    

    if protein_seq is None:
        return('Missing information in selected sequence')
    if min_pep_len > max_pep_len:
        return('Filtering not possible due to wrong lenght values')
    else:
        peptides = re.split(cleave_pattern, protein_seq)

        filtered_peptides = [pep for pep in peptides if min_pep_len <= len(pep) <= max_pep_len]

        return(filtered_peptides)

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    digest a hole protein collection

    Parameters
    protein_map: dic
    cleave_pattern: str
    min_pep_len: int
    max_pep_len: int

    Returns
    dictionary which contains list of peptides
    """

    
    protein = protein_map.keys()
    digested_collection = {}

    for pro in protein:
        digested_peptides = digest_protein_sequence(protein_map[pro], cleave_pattern, min_pep_len, max_pep_len)
        #digested_collection = ''.join(digested_peptides)
        digested_collection[pro] = digested_peptides

    return digested_collection

def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the sequence coverage (%) of a protein based on detected peptides.

    Parameters
    ----------
    protein_seq : str
        Complete protein sequence before digestion.
    peptides : list of str
        Peptide fragments detected after digestion.

    Returns
    -------
    float
        Protein coverage percentage (0–100).
    """

    # Handle empty input
    if not protein_seq or not peptides:
        return 0.0

    protein_seq = protein_seq.upper()
    coverage = set()  # store covered indices (no duplicates)

    for pep in peptides:
        pep = pep.upper()
        idx_start = 0

        # Find all occurrences of the peptide inside the protein
        while True:
            idx = protein_seq.find(pep, idx_start)
            if idx == -1:
                break

            # Add all residue indices covered by this peptide
            coverage.update(range(idx, idx + len(pep)))

            # Continue searching after this position
            idx_start = idx + 1

    # Calculate % coverage
    coverage_percent = (len(coverage) / len(protein_seq)) * 100
    return coverage_percent

def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the sequence coverage (%) of a protein based on detected peptides.

    Parameters
    ----------
    protein_seq : str
        Complete protein sequence before digestion.
    peptides : list of str
        Peptide fragments detected after digestion.

    Returns
    -------
    float
        Protein coverage percentage (0–100).
    """

    # Handle empty input
    if not protein_seq or not peptides:
        return 0.0

    protein_seq = protein_seq.upper()
    coverage = set()  # store covered indices (no duplicates)

    for pep in peptides:
        pep = pep.upper()
        idx_start = 0

        # Find all occurrences of the peptide inside the protein
        while True:
            idx = protein_seq.find(pep, idx_start)
            if idx == -1:
                break

            # Add all residue indices covered by this peptide
            coverage.update(range(idx, idx + len(pep)))

            # Continue searching after this position
            idx_start = idx + 1

    # Calculate % coverage
    coverage_percent = (len(coverage) / len(protein_seq)) * 100
    return coverage_percent