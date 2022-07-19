# Copied from https://github.com/Lattice-Automation/seqfold/blob/master/seqfold/rna.py

"""RNA enthalpy and entropy change parameters."""

from .types import Comp, MultiBranch, BpEnergy, LoopEnergy, Energies

RNA_COMPLEMENT: Comp = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}

RNA_MULTIBRANCH: MultiBranch = (2.5, 0.1, 0.4, 2.0)

RNA_NN: BpEnergy = {
    "AA/UU": (-6.8, -19.0),
    "AC/UG": (-11.4, -29.7),
    "AG/UC": (-10.5, -27.1),
    "AU/UA": (-9.4, -26.8),
    "CA/GU": (-10.4, -26.8),
    "CC/GG": (-13.4, -32.6),
    "CG/GC": (-10.6, -26.4),
    "CU/GA": (-10.5, -27.1),
    "GA/CU": (-12.4, -32.2),
    "GC/CG": (-14.9, -37.1),
    "GG/CC": (-13.4, -32.6),
    "GU/CA": (-11.4, -29.7),
    "UA/AU": (-7.7, -20.6),
    "UC/AG": (-12.4, -32.2),
    "UG/AC": (-10.4, -26.8),
    "UU/AA": (-6.8, -19.0),
}

RNA_INTERNAL_MM: BpEnergy = {
    "AA/AA": (0.0, 0.0),
    "AA/AC": (0.0, 0.0),
    "AA/AG": (0.0, 0.0),
    "AA/AU": (0.0, 0.0),
    "AA/CA": (0.0, 0.0),
    "AA/CC": (0.0, 0.0),
    "AA/CG": (0.0, 0.0),
    "AA/CU": (0.0, 0.0),
    "AA/GA": (0.0, 0.0),
    "AA/GC": (0.0, 0.0),
    "AA/GG": (0.0, 0.0),
    "AA/GU": (0.0, 0.0),
    "AA/UA": (0.0, 0.0),
    "AA/UC": (0.0, 0.0),
    "AA/UG": (0.0, 0.0),
    "AC/AA": (0.0, 0.0),
    "AC/AC": (0.0, 0.0),
    "AC/AG": (0.0, 0.0),
    "AC/AU": (0.0, 0.0),
    "AC/CA": (0.0, 0.0),
    "AC/CC": (0.0, 0.0),
    "AC/CG": (0.0, 0.0),
    "AC/CU": (0.0, 0.0),
    "AC/GA": (0.0, 0.0),
    "AC/GC": (0.0, 0.0),
    "AC/GG": (0.0, 0.0),
    "AC/GU": (0.0, 0.0),
    "AC/UA": (0.0, 0.0),
    "AC/UC": (0.0, 0.0),
    "AC/UU": (0.0, 0.0),
    "AG/AA": (0.0, 0.0),
    "AG/AC": (0.0, 0.0),
    "AG/AG": (0.0, 0.0),
    "AG/AU": (0.0, 0.0),
    "AG/CA": (0.0, 0.0),
    "AG/CC": (0.0, 0.0),
    "AG/CG": (0.0, 0.0),
    "AG/CU": (0.0, 0.0),
    "AG/GA": (0.0, 0.0),
    "AG/GC": (0.0, 0.0),
    "AG/GG": (0.0, 0.0),
    "AG/GU": (0.0, 0.0),
    "AG/UA": (0.0, 0.0),
    "AG/UG": (0.0, 0.0),
    "AG/UU": (-3.2, -8.4),
    "AU/AA": (0.0, 0.0),
    "AU/AC": (0.0, 0.0),
    "AU/AG": (0.0, 0.0),
    "AU/AU": (0.0, 0.0),
    "AU/CA": (0.0, 0.0),
    "AU/CC": (0.0, 0.0),
    "AU/CG": (0.0, 0.0),
    "AU/CU": (0.0, 0.0),
    "AU/GA": (0.0, 0.0),
    "AU/GC": (0.0, 0.0),
    "AU/GG": (0.0, 0.0),
    "AU/GU": (0.0, 0.0),
    "AU/UC": (0.0, 0.0),
    "AU/UG": (-8.8, -23.9),
    "AU/UU": (0.0, 0.0),
    "CA/AA": (0.0, 0.0),
    "CA/AC": (0.0, 0.0),
    "CA/AG": (0.0, 0.0),
    "CA/AU": (0.0, 0.0),
    "CA/CA": (0.0, 0.0),
    "CA/CC": (0.0, 0.0),
    "CA/CG": (0.0, 0.0),
    "CA/CU": (0.0, 0.0),
    "CA/GA": (0.0, 0.0),
    "CA/GC": (0.0, 0.0),
    "CA/GG": (0.0, 0.0),
    "CA/UA": (0.0, 0.0),
    "CA/UC": (0.0, 0.0),
    "CA/UG": (0.0, 0.0),
    "CA/UU": (0.0, 0.0),
    "CC/AA": (0.0, 0.0),
    "CC/AC": (0.0, 0.0),
    "CC/AG": (0.0, 0.0),
    "CC/AU": (0.0, 0.0),
    "CC/CA": (0.0, 0.0),
    "CC/CC": (0.0, 0.0),
    "CC/CG": (0.0, 0.0),
    "CC/CU": (0.0, 0.0),
    "CC/GA": (0.0, 0.0),
    "CC/GC": (0.0, 0.0),
    "CC/GU": (0.0, 0.0),
    "CC/UA": (0.0, 0.0),
    "CC/UC": (0.0, 0.0),
    "CC/UG": (0.0, 0.0),
    "CC/UU": (0.0, 0.0),
    "CG/AA": (0.0, 0.0),
    "CG/AC": (0.0, 0.0),
    "CG/AG": (0.0, 0.0),
    "CG/AU": (0.0, 0.0),
    "CG/CA": (0.0, 0.0),
    "CG/CC": (0.0, 0.0),
    "CG/CG": (0.0, 0.0),
    "CG/CU": (0.0, 0.0),
    "CG/GA": (0.0, 0.0),
    "CG/GG": (0.0, 0.0),
    "CG/GU": (-5.6, -13.5),
    "CG/UA": (0.0, 0.0),
    "CG/UC": (0.0, 0.0),
    "CG/UG": (0.0, 0.0),
    "CG/UU": (0.0, 0.0),
    "CU/AA": (0.0, 0.0),
    "CU/AC": (0.0, 0.0),
    "CU/AG": (0.0, 0.0),
    "CU/AU": (0.0, 0.0),
    "CU/CA": (0.0, 0.0),
    "CU/CC": (0.0, 0.0),
    "CU/CG": (0.0, 0.0),
    "CU/CU": (0.0, 0.0),
    "CU/GC": (0.0, 0.0),
    "CU/GG": (-12.1, -32.2),
    "CU/GU": (0.0, 0.0),
    "CU/UA": (0.0, 0.0),
    "CU/UC": (0.0, 0.0),
    "CU/UG": (0.0, 0.0),
    "CU/UU": (0.0, 0.0),
    "GA/AA": (0.0, 0.0),
    "GA/AC": (0.0, 0.0),
    "GA/AG": (0.0, 0.0),
    "GA/AU": (0.0, 0.0),
    "GA/CA": (0.0, 0.0),
    "GA/CC": (0.0, 0.0),
    "GA/CG": (0.0, 0.0),
    "GA/GA": (0.0, 0.0),
    "GA/GC": (0.0, 0.0),
    "GA/GG": (0.0, 0.0),
    "GA/GU": (0.0, 0.0),
    "GA/UA": (0.0, 0.0),
    "GA/UC": (0.0, 0.0),
    "GA/UG": (0.0, 0.0),
    "GA/UU": (-12.8, -37.1),
    "GC/AA": (0.0, 0.0),
    "GC/AC": (0.0, 0.0),
    "GC/AG": (0.0, 0.0),
    "GC/AU": (0.0, 0.0),
    "GC/CA": (0.0, 0.0),
    "GC/CC": (0.0, 0.0),
    "GC/CU": (0.0, 0.0),
    "GC/GA": (0.0, 0.0),
    "GC/GC": (0.0, 0.0),
    "GC/GG": (0.0, 0.0),
    "GC/GU": (0.0, 0.0),
    "GC/UA": (0.0, 0.0),
    "GC/UC": (0.0, 0.0),
    "GC/UG": (-12.6, -32.6),
    "GC/UU": (0.0, 0.0),
    "GG/AA": (0.0, 0.0),
    "GG/AC": (0.0, 0.0),
    "GG/AG": (0.0, 0.0),
    "GG/AU": (0.0, 0.0),
    "GG/CA": (0.0, 0.0),
    "GG/CG": (0.0, 0.0),
    "GG/CU": (-8.3, -21.9),
    "GG/GA": (0.0, 0.0),
    "GG/GC": (0.0, 0.0),
    "GG/GG": (0.0, 0.0),
    "GG/GU": (0.0, 0.0),
    "GG/UA": (0.0, 0.0),
    "GG/UC": (-12.1, -32.2),
    "GG/UG": (0.0, 0.0),
    "GG/UU": (-13.5, -41.9),
    "GU/AA": (0.0, 0.0),
    "GU/AC": (0.0, 0.0),
    "GU/AG": (0.0, 0.0),
    "GU/AU": (0.0, 0.0),
    "GU/CC": (0.0, 0.0),
    "GU/CG": (-12.6, -32.6),
    "GU/CU": (0.0, 0.0),
    "GU/GA": (0.0, 0.0),
    "GU/GC": (0.0, 0.0),
    "GU/GG": (0.0, 0.0),
    "GU/GU": (0.0, 0.0),
    "GU/UA": (-8.8, -23.9),
    "GU/UC": (0.0, 0.0),
    "GU/UG": (-14.6, -51.3),
    "GU/UU": (0.0, 0.0),
    "UA/AA": (0.0, 0.0),
    "UA/AC": (0.0, 0.0),
    "UA/AG": (0.0, 0.0),
    "UA/CA": (0.0, 0.0),
    "UA/CC": (0.0, 0.0),
    "UA/CG": (0.0, 0.0),
    "UA/CU": (0.0, 0.0),
    "UA/GA": (0.0, 0.0),
    "UA/GC": (0.0, 0.0),
    "UA/GG": (0.0, 0.0),
    "UA/GU": (-7.0, -19.3),
    "UA/UA": (0.0, 0.0),
    "UA/UC": (0.0, 0.0),
    "UA/UG": (0.0, 0.0),
    "UA/UU": (0.0, 0.0),
    "UC/AA": (0.0, 0.0),
    "UC/AC": (0.0, 0.0),
    "UC/AU": (0.0, 0.0),
    "UC/CA": (0.0, 0.0),
    "UC/CC": (0.0, 0.0),
    "UC/CG": (0.0, 0.0),
    "UC/CU": (0.0, 0.0),
    "UC/GA": (0.0, 0.0),
    "UC/GC": (0.0, 0.0),
    "UC/GG": (-8.3, -21.9),
    "UC/GU": (0.0, 0.0),
    "UC/UA": (0.0, 0.0),
    "UC/UC": (0.0, 0.0),
    "UC/UG": (0.0, 0.0),
    "UC/UU": (0.0, 0.0),
    "UG/AA": (0.0, 0.0),
    "UG/AG": (0.0, 0.0),
    "UG/AU": (-7.0, -19.3),
    "UG/CA": (0.0, 0.0),
    "UG/CC": (0.0, 0.0),
    "UG/CG": (0.0, 0.0),
    "UG/CU": (0.0, 0.0),
    "UG/GA": (0.0, 0.0),
    "UG/GC": (-5.6, -13.5),
    "UG/GG": (0.0, 0.0),
    "UG/GU": (-9.3, -31.0),
    "UG/UA": (0.0, 0.0),
    "UG/UC": (0.0, 0.0),
    "UG/UG": (0.0, 0.0),
    "UG/UU": (0.0, 0.0),
    "UU/AC": (0.0, 0.0),
    "UU/AG": (-12.8, -37.1),
    "UU/AU": (0.0, 0.0),
    "UU/CA": (0.0, 0.0),
    "UU/CC": (0.0, 0.0),
    "UU/CG": (0.0, 0.0),
    "UU/CU": (0.0, 0.0),
    "UU/GA": (-3.2, -8.4),
    "UU/GC": (0.0, 0.0),
    "UU/GG": (-13.5, -41.9),
    "UU/GU": (0.0, 0.0),
    "UU/UA": (0.0, 0.0),
    "UU/UC": (0.0, 0.0),
    "UU/UG": (0.0, 0.0),
    "UU/UU": (0.0, 0.0),
}

RNA_TERMINAL_MM: BpEnergy = {
    "AA/AA": (0.0, 0.0),
    "AA/AC": (0.0, 0.0),
    "AA/AG": (0.0, 0.0),
    "AA/AU": (0.0, 0.0),
    "AA/CA": (0.0, 0.0),
    "AA/CC": (0.0, 0.0),
    "AA/CG": (0.0, 0.0),
    "AA/CU": (0.0, 0.0),
    "AA/GA": (0.0, 0.0),
    "AA/GC": (0.0, 0.0),
    "AA/GG": (0.0, 0.0),
    "AA/GU": (0.0, 0.0),
    "AA/UA": (-3.9, -10.0),
    "AA/UC": (2.0, 9.7),
    "AA/UG": (-3.5, -8.7),
    "AA/UU": (2.0, 9.7),
    "AC/AA": (0.0, 0.0),
    "AC/AC": (0.0, 0.0),
    "AC/AG": (0.0, 0.0),
    "AC/AU": (0.0, 0.0),
    "AC/CA": (0.0, 0.0),
    "AC/CC": (0.0, 0.0),
    "AC/CG": (0.0, 0.0),
    "AC/CU": (0.0, 0.0),
    "AC/GA": (0.0, 0.0),
    "AC/GC": (0.0, 0.0),
    "AC/GG": (0.0, 0.0),
    "AC/GU": (0.0, 0.0),
    "AC/UA": (-2.3, -5.5),
    "AC/UC": (6.0, 21.6),
    "AC/UG": (-2.3, -5.5),
    "AC/UU": (-0.3, 1.3),
    "AG/AA": (0.0, 0.0),
    "AG/AC": (0.0, 0.0),
    "AG/AG": (0.0, 0.0),
    "AG/AU": (0.0, 0.0),
    "AG/CA": (0.0, 0.0),
    "AG/CC": (0.0, 0.0),
    "AG/CG": (0.0, 0.0),
    "AG/CU": (0.0, 0.0),
    "AG/GA": (0.0, 0.0),
    "AG/GC": (0.0, 0.0),
    "AG/GG": (0.0, 0.0),
    "AG/GU": (0.0, 0.0),
    "AG/UA": (-3.1, -7.4),
    "AG/UC": (2.0, 9.7),
    "AG/UG": (-3.5, -8.7),
    "AG/UU": (2.0, 9.7),
    "AU/AA": (0.0, 0.0),
    "AU/AC": (0.0, 0.0),
    "AU/AG": (0.0, 0.0),
    "AU/AU": (0.0, 0.0),
    "AU/CA": (0.0, 0.0),
    "AU/CC": (0.0, 0.0),
    "AU/CG": (0.0, 0.0),
    "AU/CU": (0.0, 0.0),
    "AU/GA": (0.0, 0.0),
    "AU/GC": (0.0, 0.0),
    "AU/GG": (0.0, 0.0),
    "AU/GU": (0.0, 0.0),
    "AU/UA": (-2.3, -5.5),
    "AU/UC": (4.6, 17.4),
    "AU/UG": (-2.3, -5.5),
    "AU/UU": (-1.7, -2.9),
    "CA/AA": (0.0, 0.0),
    "CA/AC": (0.0, 0.0),
    "CA/AG": (0.0, 0.0),
    "CA/AU": (0.0, 0.0),
    "CA/CA": (0.0, 0.0),
    "CA/CC": (0.0, 0.0),
    "CA/CG": (0.0, 0.0),
    "CA/CU": (0.0, 0.0),
    "CA/GA": (-9.1, -24.5),
    "CA/GC": (-5.6, -13.2),
    "CA/GG": (-5.6, -13.5),
    "CA/GU": (-5.6, -13.2),
    "CA/UA": (0.0, 0.0),
    "CA/UC": (0.0, 0.0),
    "CA/UG": (0.0, 0.0),
    "CA/UU": (0.0, 0.0),
    "CC/AA": (0.0, 0.0),
    "CC/AC": (0.0, 0.0),
    "CC/AG": (0.0, 0.0),
    "CC/AU": (0.0, 0.0),
    "CC/CA": (0.0, 0.0),
    "CC/CC": (0.0, 0.0),
    "CC/CG": (0.0, 0.0),
    "CC/CU": (0.0, 0.0),
    "CC/GA": (-5.7, -15.2),
    "CC/GC": (-3.4, -7.4),
    "CC/GG": (-5.7, -15.2),
    "CC/GU": (-2.7, -6.1),
    "CC/UA": (0.0, 0.0),
    "CC/UC": (0.0, 0.0),
    "CC/UG": (0.0, 0.0),
    "CC/UU": (0.0, 0.0),
    "CG/AA": (0.0, 0.0),
    "CG/AC": (0.0, 0.0),
    "CG/AG": (0.0, 0.0),
    "CG/AU": (0.0, 0.0),
    "CG/CA": (0.0, 0.0),
    "CG/CC": (0.0, 0.0),
    "CG/CG": (0.0, 0.0),
    "CG/CU": (0.0, 0.0),
    "CG/GA": (-8.2, -21.9),
    "CG/GC": (-5.6, -13.2),
    "CG/GG": (-9.2, -24.5),
    "CG/GU": (-5.6, -13.2),
    "CG/UA": (0.0, 0.0),
    "CG/UC": (0.0, 0.0),
    "CG/UG": (0.0, 0.0),
    "CG/UU": (0.0, 0.0),
    "CU/AA": (0.0, 0.0),
    "CU/AC": (0.0, 0.0),
    "CU/AG": (0.0, 0.0),
    "CU/AU": (0.0, 0.0),
    "CU/CA": (0.0, 0.0),
    "CU/CC": (0.0, 0.0),
    "CU/CG": (0.0, 0.0),
    "CU/CU": (0.0, 0.0),
    "CU/GA": (-5.7, -15.2),
    "CU/GC": (-5.3, -12.6),
    "CU/GG": (-5.7, -15.2),
    "CU/GU": (-8.6, -23.9),
    "CU/UA": (0.0, 0.0),
    "CU/UC": (0.0, 0.0),
    "CU/UG": (0.0, 0.0),
    "CU/UU": (0.0, 0.0),
    "GA/AA": (0.0, 0.0),
    "GA/AC": (0.0, 0.0),
    "GA/AG": (0.0, 0.0),
    "GA/AU": (0.0, 0.0),
    "GA/CA": (-5.2, -13.2),
    "GA/CC": (-4.0, -8.1),
    "GA/CG": (-5.6, -13.9),
    "GA/CU": (-4.0, -8.1),
    "GA/GA": (0.0, 0.0),
    "GA/GC": (0.0, 0.0),
    "GA/GG": (0.0, 0.0),
    "GA/GU": (0.0, 0.0),
    "GA/UA": (-3.4, -10.0),
    "GA/UC": (2.0, 9.7),
    "GA/UG": (-3.5, -8.7),
    "GA/UU": (2.0, 9.7),
    "GC/AA": (0.0, 0.0),
    "GC/AC": (0.0, 0.0),
    "GC/AG": (0.0, 0.0),
    "GC/AU": (0.0, 0.0),
    "GC/CA": (-7.2, -19.7),
    "GC/CC": (0.5, 3.9),
    "GC/CG": (-7.2, -19.7),
    "GC/CU": (-4.2, -11.9),
    "GC/GA": (0.0, 0.0),
    "GC/GC": (0.0, 0.0),
    "GC/GG": (0.0, 0.0),
    "GC/GU": (0.0, 0.0),
    "GC/UA": (-2.3, -5.5),
    "GC/UC": (6.0, 21.6),
    "GC/UG": (-2.3, -5.5),
    "GC/UU": (-0.3, 1.3),
    "GG/AA": (0.0, 0.0),
    "GG/AC": (0.0, 0.0),
    "GG/AG": (0.0, 0.0),
    "GG/AU": (0.0, 0.0),
    "GG/CA": (-7.1, -17.7),
    "GG/CC": (-4.0, -8.1),
    "GG/CG": (-6.2, -15.5),
    "GG/CU": (-4.0, -8.1),
    "GG/GA": (0.0, 0.0),
    "GG/GC": (0.0, 0.0),
    "GG/GG": (0.0, 0.0),
    "GG/GU": (0.0, 0.0),
    "GG/UA": (-0.6, 0.0),
    "GG/UC": (2.0, 9.7),
    "GG/UG": (-3.5, -8.7),
    "GG/UU": (2.0, 9.7),
    "GU/AA": (0.0, 0.0),
    "GU/AC": (0.0, 0.0),
    "GU/AG": (0.0, 0.0),
    "GU/AU": (0.0, 0.0),
    "GU/CA": (-7.2, -19.7),
    "GU/CC": (-0.3, 2.3),
    "GU/CG": (-7.2, -19.7),
    "GU/CU": (-5.0, -13.9),
    "GU/GA": (0.0, 0.0),
    "GU/GC": (0.0, 0.0),
    "GU/GG": (0.0, 0.0),
    "GU/GU": (0.0, 0.0),
    "GU/UA": (-2.3, -5.5),
    "GU/UC": (4.6, 17.4),
    "GU/UG": (-2.3, -5.5),
    "GU/UU": (1.6, 7.1),
    "UA/AA": (-4.0, -9.7),
    "UA/AC": (-6.3, -17.7),
    "UA/AG": (-8.9, -25.1),
    "UA/AU": (-6.3, -17.7),
    "UA/CA": (0.0, 0.0),
    "UA/CC": (0.0, 0.0),
    "UA/CG": (0.0, 0.0),
    "UA/CU": (0.0, 0.0),
    "UA/GA": (-4.8, -12.3),
    "UA/GC": (-6.3, -17.7),
    "UA/GG": (-8.9, -25.1),
    "UA/GU": (-6.3, -17.7),
    "UA/UA": (0.0, 0.0),
    "UA/UC": (0.0, 0.0),
    "UA/UG": (0.0, 0.0),
    "UA/UU": (0.0, 0.0),
    "UC/AA": (-4.3, -11.6),
    "UC/AC": (-5.1, -14.5),
    "UC/AG": (-4.3, -11.6),
    "UC/AU": (-1.8, -4.2),
    "UC/CA": (0.0, 0.0),
    "UC/CC": (0.0, 0.0),
    "UC/CG": (0.0, 0.0),
    "UC/CU": (0.0, 0.0),
    "UC/GA": (-4.3, -11.6),
    "UC/GC": (-5.1, -14.5),
    "UC/GG": (-4.3, -11.6),
    "UC/GU": (-1.8, -4.2),
    "UC/UA": (0.0, 0.0),
    "UC/UC": (0.0, 0.0),
    "UC/UG": (0.0, 0.0),
    "UC/UU": (0.0, 0.0),
    "UG/AA": (-3.8, -8.7),
    "UG/AC": (-6.3, -17.7),
    "UG/AG": (-8.9, -24.8),
    "UG/AU": (-6.3, -17.7),
    "UG/CA": (0.0, 0.0),
    "UG/CC": (0.0, 0.0),
    "UG/CG": (0.0, 0.0),
    "UG/CU": (0.0, 0.0),
    "UG/GA": (3.1, 11.6),
    "UG/GC": (-6.3, -17.7),
    "UG/GG": (-1.5, -2.3),
    "UG/GU": (-6.3, -17.7),
    "UG/UA": (0.0, 0.0),
    "UG/UC": (0.0, 0.0),
    "UG/UG": (0.0, 0.0),
    "UG/UU": (0.0, 0.0),
    "UU/AA": (-4.3, -11.6),
    "UU/AC": (-1.4, -2.6),
    "UU/AG": (-4.3, -11.6),
    "UU/AU": (1.4, 6.1),
    "UU/CA": (0.0, 0.0),
    "UU/CC": (0.0, 0.0),
    "UU/CG": (0.0, 0.0),
    "UU/CU": (0.0, 0.0),
    "UU/GA": (-4.3, -11.6),
    "UU/GC": (-1.4, -2.6),
    "UU/GG": (-4.3, -11.6),
    "UU/GU": (1.4, 6.1),
    "UU/UA": (0.0, 0.0),
    "UU/UC": (0.0, 0.0),
    "UU/UG": (0.0, 0.0),
    "UU/UU": (0.0, 0.0),
}

RNA_DE: BpEnergy = {
    "AA/A.": (0.0, 0.0),
    "AC/A.": (0.0, 0.0),
    "AG/A.": (0.0, 0.0),
    "AU/A.": (0.0, 0.0),
    "AA/C.": (0.0, 0.0),
    "AC/C.": (0.0, 0.0),
    "AG/C.": (0.0, 0.0),
    "AU/C.": (0.0, 0.0),
    "AA/G.": (0.0, 0.0),
    "AC/G.": (0.0, 0.0),
    "AG/G.": (0.0, 0.0),
    "AU/G.": (0.0, 0.0),
    "AA/U.": (-4.9, -13.2),
    "AC/U.": (-0.9, -1.3),
    "AG/U.": (-5.5, -15.2),
    "AU/U.": (-2.3, -5.5),
    "CA/A.": (0.0, 0.0),
    "CC/A.": (0.0, 0.0),
    "CG/A.": (0.0, 0.0),
    "CU/A.": (0.0, 0.0),
    "CA/C.": (0.0, 0.0),
    "CC/C.": (0.0, 0.0),
    "CG/C.": (0.0, 0.0),
    "CU/C.": (0.0, 0.0),
    "CA/G.": (-9.0, -23.5),
    "CC/G.": (-4.1, -10.6),
    "CG/G.": (-8.6, -22.2),
    "CU/G.": (-7.5, -20.3),
    "CA/U.": (0.0, 0.0),
    "CC/U.": (0.0, 0.0),
    "CG/U.": (0.0, 0.0),
    "CU/U.": (0.0, 0.0),
    "GA/A.": (0.0, 0.0),
    "GC/A.": (0.0, 0.0),
    "GG/A.": (0.0, 0.0),
    "GU/A.": (0.0, 0.0),
    "GA/C.": (-7.4, -20.3),
    "GC/C.": (-2.8, -7.7),
    "GG/C.": (-6.4, -16.4),
    "GU/C.": (-3.6, -9.7),
    "GA/G.": (0.0, 0.0),
    "GC/G.": (0.0, 0.0),
    "GG/G.": (0.0, 0.0),
    "GU/G.": (0.0, 0.0),
    "GA/U.": (-4.9, -13.2),
    "GC/U.": (-0.9, -1.3),
    "GG/U.": (-5.5, -15.2),
    "GU/U.": (-2.3, -5.5),
    "UA/A.": (-5.7, -16.1),
    "UC/A.": (-0.7, -1.9),
    "UG/A.": (-5.8, -16.4),
    "UU/A.": (-2.2, -6.8),
    "UA/C.": (0.0, 0.0),
    "UC/C.": (0.0, 0.0),
    "UG/C.": (0.0, 0.0),
    "UU/C.": (0.0, 0.0),
    "UA/G.": (-5.7, -16.1),
    "UC/G.": (-0.7, -1.9),
    "UG/G.": (-5.8, -16.4),
    "UU/G.": (-2.2, -6.8),
    "UA/U.": (0.0, 0.0),
    "UC/U.": (0.0, 0.0),
    "UG/U.": (0.0, 0.0),
    "UU/U.": (0.0, 0.0),
    "A./AA": (0.0, 0.0),
    "A./AC": (0.0, 0.0),
    "A./AG": (0.0, 0.0),
    "A./AU": (0.0, 0.0),
    "A./CA": (0.0, 0.0),
    "A./CC": (0.0, 0.0),
    "A./CG": (0.0, 0.0),
    "A./CU": (0.0, 0.0),
    "A./GA": (0.0, 0.0),
    "A./GC": (0.0, 0.0),
    "A./GG": (0.0, 0.0),
    "A./GU": (0.0, 0.0),
    "A./UA": (-0.5, -0.6),
    "A./UC": (6.9, 22.6),
    "A./UG": (0.6, 2.6),
    "A./UU": (0.6, 2.6),
    "C./AA": (0.0, 0.0),
    "C./AC": (0.0, 0.0),
    "C./AG": (0.0, 0.0),
    "C./AU": (0.0, 0.0),
    "C./CA": (0.0, 0.0),
    "C./CC": (0.0, 0.0),
    "C./CG": (0.0, 0.0),
    "C./CU": (0.0, 0.0),
    "C./GA": (-1.6, -4.5),
    "C./GC": (0.7, 3.2),
    "C./GG": (-4.6, -14.8),
    "C./GU": (-0.4, -1.3),
    "C./UA": (0.0, 0.0),
    "C./UC": (0.0, 0.0),
    "C./UG": (0.0, 0.0),
    "C./UU": (0.0, 0.0),
    "G./AA": (0.0, 0.0),
    "G./AC": (0.0, 0.0),
    "G./AG": (0.0, 0.0),
    "G./AU": (0.0, 0.0),
    "G./CA": (-2.4, -6.1),
    "G./CC": (3.3, 11.6),
    "G./CG": (0.8, 3.2),
    "G./CU": (-1.4, -4.2),
    "G./GA": (0.0, 0.0),
    "G./GC": (0.0, 0.0),
    "G./GG": (0.0, 0.0),
    "G./GU": (0.0, 0.0),
    "G./UA": (-0.5, -0.6),
    "G./UC": (6.9, 22.6),
    "G./UG": (0.6, 2.6),
    "G./UU": (0.6, 2.6),
    "U./AA": (1.6, 6.1),
    "U./AC": (2.2, 8.1),
    "U./AG": (0.7, 3.5),
    "U./AU": (3.1, 10.6),
    "U./CA": (0.0, 0.0),
    "U./CC": (0.0, 0.0),
    "U./CG": (0.0, 0.0),
    "U./CU": (0.0, 0.0),
    "U./GA": (1.6, 6.1),
    "U./GC": (2.2, 8.1),
    "U./GG": (0.7, 3.5),
    "U./GU": (3.1, 10.6),
    "U./UA": (0.0, 0.0),
    "U./UC": (0.0, 0.0),
    "U./UG": (0.0, 0.0),
    "U./UU": (0.0, 0.0),
}

RNA_INTERNAL_LOOPS: LoopEnergy = {
    1: (0.0, 0.0),
    2: (0.0, 0.0),
    3: (0.0, 0.0),
    4: (-7.2, -26.8),
    5: (-6.8, -28.4),
    6: (-1.3, -10.6),
    7: (-1.3, -11.0),
    8: (-1.3, -11.6),
    9: (-1.3, -11.9),
    10: (-1.3, -12.3),
    11: (-1.3, -12.6),
    12: (-1.3, -12.9),
    13: (-1.3, -13.2),
    14: (-1.3, -13.5),
    15: (-1.3, -13.5),
    16: (-1.3, -13.9),
    17: (-1.3, -14.2),
    18: (-1.3, -14.2),
    19: (-1.3, -14.5),
    20: (-1.3, -14.8),
    21: (-1.3, -14.8),
    22: (-1.3, -15.2),
    23: (-1.3, -15.2),
    24: (-1.3, -15.5),
    25: (-1.3, -15.5),
    26: (-1.3, -15.5),
    27: (-1.3, -15.8),
    28: (-1.3, -15.8),
    29: (-1.3, -16.1),
    30: (-1.3, -16.1),
}

RNA_BULGE_LOOPS: LoopEnergy = {
    1: (10.6, 21.9),
    2: (7.1, 13.9),
    3: (7.1, 12.6),
    4: (7.1, 11.3),
    5: (7.1, 10.0),
    6: (7.1, 8.7),
    7: (7.1, 8.1),
    8: (7.1, 7.7),
    9: (7.1, 7.4),
    10: (7.1, 7.1),
    11: (7.1, 6.8),
    12: (7.1, 6.4),
    13: (7.1, 6.1),
    14: (7.1, 5.8),
    15: (7.1, 5.5),
    16: (7.1, 5.5),
    17: (7.1, 5.2),
    18: (7.1, 5.2),
    19: (7.1, 4.8),
    20: (7.1, 4.5),
    21: (7.1, 4.5),
    22: (7.1, 4.2),
    23: (7.1, 4.2),
    24: (7.1, 4.2),
    25: (7.1, 3.9),
    26: (7.1, 3.9),
    27: (7.1, 3.5),
    28: (7.1, 3.5),
    29: (7.1, 3.5),
    30: (7.1, 3.2),
}

RNA_HAIRPIN_LOOPS: LoopEnergy = {
    1: (0.0, 0.0),
    2: (0.0, 0.0),
    3: (1.3, -13.2),
    4: (4.8, -2.6),
    5: (3.6, -6.8),
    6: (-2.9, -26.8),
    7: (1.3, -15.2),
    8: (-2.9, -27.1),
    9: (5.0, -4.5),
    10: (5.0, -4.8),
    11: (5.0, -5.2),
    12: (5.0, -5.5),
    13: (5.0, -5.8),
    14: (5.0, -6.1),
    15: (5.0, -6.1),
    16: (5.0, -6.4),
    17: (5.0, -6.8),
    18: (5.0, -6.8),
    19: (5.0, -7.1),
    20: (5.0, -7.1),
    21: (5.0, -7.4),
    22: (5.0, -7.4),
    23: (5.0, -7.7),
    24: (5.0, -7.7),
    25: (5.0, -8.1),
    26: (5.0, -8.1),
    27: (5.0, -8.1),
    28: (5.0, -8.4),
    29: (5.0, -8.4),
    30: (5.0, -8.7),
}


# the energies are the same for each loop stack in the
# reverse complementary direction
RNA_NN.update({k[::-1]: v for k, v in RNA_NN.items()})
RNA_INTERNAL_MM.update(
    {k[::-1]: v for k, v in RNA_INTERNAL_MM.items() if k[::-1] not in RNA_INTERNAL_MM}
)
RNA_TERMINAL_MM.update(
    {k[::-1]: v for k, v in RNA_TERMINAL_MM.items() if k[::-1] not in RNA_TERMINAL_MM}
)
RNA_DE.update({k[::-1]: v for k, v in RNA_DE.items() if k[::-1] not in RNA_DE})

RNA_ENERGIES = Energies(
    RNA_BULGE_LOOPS,
    RNA_COMPLEMENT,
    RNA_DE,
    RNA_HAIRPIN_LOOPS,
    RNA_MULTIBRANCH,
    RNA_INTERNAL_LOOPS,
    RNA_INTERNAL_MM,
    RNA_NN,
    RNA_TERMINAL_MM,
    None,
)
