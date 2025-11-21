from Proteosim_EA.liquid_chromatography import predict_lc_retention_times
from Proteosim_EA.liquid_chromatography import select_retention_time_window
from Proteosim_EA.liquid_chromatography import achrom

def test_predict_lc_retention_times():
    test_peptides = ["AAAAAA", "PPPPPP", "CCCCCC"]

    expected = {} #dictionary, das Peptidsequenzen als key und vorhergesagten RT als value enthält
    for pep in test_peptides:
        expected[pep] = round(achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0), 2) #ausgegebene RT auf zwei Nachkommastellen gerundet

    # zu testende Funktion
    actual = predict_lc_retention_times(test_peptides)

    assert actual == expected


def test_test_select_retention_time_window():

    test_peptide_rt_map = {
        "PEP1": 12.3,
        "PEP2": 25.7,
        "PEP3": 8.9,
        "PEP4": 18.5
    }
    
    selected = select_retention_time_window(
        test_peptide_rt_map,
        lower_ret_time= 10,
        upper_ret_time= 20
    )

    expected = {"PEP1": 12.3, "PEP4":18.5}

    assert selected == expected