from Proteosim_EA.mass_spectra_simulation import calculate_mol_mass
from Proteosim_EA.mass_spectra_simulation import calculate_mol_mass_collection
from Proteosim_EA.mass_spectra_simulation import calculate_mz_collection
from Proteosim_EA.mass_spectra_simulation import fragment_peptide
from Proteosim_EA.mass_spectra_simulation import amino_acid_mass_dalton

def test_calculate_mol_mass():
    aa_mass_dict = {"A": 10, "C": 20, "D": 30, "C": 40, "E": 50} 
    test_peptide = "ACDE"

    
    expected_mass = {"ACDE": aa_mass_dict["A"] + aa_mass_dict["C"] + aa_mass_dict["D"] + aa_mass_dict["E"]}

    
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)

    
    assert calculated_mass == expected_mass

def test_calculate_mol_mass_collection():
    
    aa_mass_dict = {"A": 10, "C": 20, "D": 30, "C": 40, "E": 50}

    test_peptides = ["ACDE", "AAA"]

    expected = {
        "ACDE": aa_mass_dict["A"] + aa_mass_dict["C"] + aa_mass_dict["D"] + aa_mass_dict["E"],
        "AAA": aa_mass_dict["A"] * 3
    }

    actual = calculate_mol_mass_collection(test_peptides, aa_mass_dict)

    assert actual == expected

def test_calculate_mz_collection():
    mass_ACDE = list(calculate_mol_mass('ACDE', amino_acid_mass_dalton).values())[0]
    mass_AAA = list(calculate_mol_mass('AAA', amino_acid_mass_dalton).values())[0]

    
    test_peptide_mass_map = {
        "ACDE": mass_ACDE,
        "AAA": mass_AAA
    }

    expected = {
        "ACDE": (mass_ACDE + 2*1.007)/2,
        "AAA": (mass_AAA + 2*1.007)/2
    }

    
    actual = calculate_mz_collection(test_peptide_mass_map, charge=2)

    
    assert actual == expected

def test_fragment_peptide():
    peptide = 'PEPT'
    expected = ['P', 'PE', 'PEP', 'T', 'PT', 'EPT']
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)

