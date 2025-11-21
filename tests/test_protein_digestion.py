from Proteosim_EA.protein_digestion import digest_protein_collection
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

