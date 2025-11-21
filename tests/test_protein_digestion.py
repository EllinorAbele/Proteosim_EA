from Proteosim_EA.protein_digestion import digest_protein_collection
from Proteosim_EA.protein_digestion import compute_sequence_coverage
def test_test_digest_protein_collection():
    dummy_proteins = {
        "protein1": "APCDERFGFFFHACDRPEFGH",
        #"protein2": "KKKPKKKP"
    }
    
    enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)', # nach Lysin geschnitten
    'LysN': r'(?=K)',  # vor Lysin geschnitten
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',       # eckige Klammer ist entweder oder
    # schneidet am C-Terminus von Lysin ODER Arginin, und Befehl noch weiter "vor" zu schauen, nur wenn da kein Prolin drauf folgt
    }

    test_digested_peptides_collection = digest_protein_collection(
        dummy_proteins,
        cleave_pattern=enzyme_cleavage_patterns["Trypsin"],
        min_pep_len=5,
        max_pep_len=30
    )

    assert test_digested_peptides_collection.get("protein1") == ["APCDER", "FGFFFHACDRPEFGH"]
    # assert test_digested_peptides_collection['protein2'] == [...]

def test_test_compute_sequence_coverage():
    # Case 1 — 0% coverage (no peptide matches)
    dummy_prot_seq = "ABCDEFGH"
    dummy_peps = ["XYZ", "LMN"]
    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 0.0

    # Case 2 — 100% coverage (full match with overlapping peptides)
    dummy_prot_seq = "ABCDEFGH"
    dummy_peps = ["ABC", "DEF", "GH"]
    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 100.0

    # Case 3 — partial coverage (~50%)
    dummy_prot_seq = "ABCDEFGH"
    dummy_peps = ["ABC", "GH"]  # covers 3 + 2 = 5 positions out of 8
    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == (5 / 8) * 100


