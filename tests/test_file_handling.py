from Proteosim_EA.file_handling import read_fasta

def test_test_fasta_reader():
    tmp_fasta_path = 'data/dummy_proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary

    # überprüft, ob die Proteine in der Testsequenz in der Dummy File 
    assert protein_map["PROTEIN_ID_1"] == "DUMMYSEQUENCEONE"
    assert protein_map["PROTEIN_ID_2"] == "DUMMYSEQUENCETWO"
