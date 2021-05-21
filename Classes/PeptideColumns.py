from typing import Tuple
from .BaseClasses.ColumnNames import ColumnNames, dataclass


@dataclass
class PeptideColumns(ColumnNames):

    """Хранит номера и имена столбцов для Peptide таблиц

    Attributes:
        accession: номер и имя столбца Accession
        sc: номер и имя столбца Sc
        precursorSignal: номер и имя столбца Precursor Signal
        sequence: номер и имя столбца Sequence
        confidence: номер и имя столбца Confidence
    """

    accession: Tuple[int, str] = (3, "Accessions")
    sc: Tuple[int, str] = (20, "Score")
    precursorSignal: Tuple[int, str] = (23, "Intensity (Peptide)")
    sequence: Tuple[int, str] = (8, "Sequence")
    confidence: Tuple[int, str] = (6, "Best Conf (Peptide)")
