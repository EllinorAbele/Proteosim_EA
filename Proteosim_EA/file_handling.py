def read_fasta(filepath):
    """
    Reads a FASTA file and returns a dictionary mapping protein IDs to their sequences.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Dictionary where keys are protein IDs and values are the full amino-acid sequences as strings.
    
    """
    protein_map = {}
    current_id = None
    current_sequence = []

    with open(filepath, 'r', encoding='utf-8') as fasta_handle:
        for line in fasta_handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('>'):
                if current_id is not None:
                    protein_map[current_id] = ''.join(current_sequence)
                    current_sequence = []
                current_id = stripped.split('|')[1]
            else:
                current_sequence.append(stripped)

    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)

    return protein_map

