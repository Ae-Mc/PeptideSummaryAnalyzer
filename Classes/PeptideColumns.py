from typing import Tuple
from .BaseClasses.ColumnNames import ColumnNames, dataclass


@dataclass
class PeptideColumns(ColumnNames):

    """Хранит номера и имена столбцов для Peptide таблиц

    Attributes:
        accession: номер и имя столбца Accession
        unused: номер и имя столбца Unused
        sc: номер и имя столбца Sc
        precursorSignal: номер и имя столбца Precursor Signal
        sequence: номер и имя столбца Sequence
        confidence: номер и имя столбца Confidence
        contribution: номер и имя столбца Contribution
    """

    accession: Tuple[int, str] = (6, "Accessions")
    unused: Tuple[int, str] = (1, "Unused")
    sc: Tuple[int, str] = (21, "Sc")
    precursorSignal: Tuple[int, str] = (24, "PrecursorSignal")
    sequence: Tuple[int, str] = (12, "Sequence")
    confidence: Tuple[int, str] = (11, "Conf")
    contribution: Tuple[int, str] = (10, "Contrib")
