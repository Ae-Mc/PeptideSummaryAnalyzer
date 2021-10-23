"""См. класс PeptideColumns"""
from typing import Tuple
from classes.base_classes.column_names import ColumnNames, dataclass


@dataclass
class PeptideColumns(ColumnNames):

    """Хранит номера и имена столбцов для Peptide таблиц

    Attributes:
        accession: номер и имя столбца Accession
        score: номер и имя столбца Score
        peptide_intensity: номер и имя столбца Peptide Intensity
        sequence: номер и имя столбца Sequence
        confidence: номер и имя столбца Confidence
    """

    accession: Tuple[int, str] = (3, "Accessions")
    score: Tuple[int, str] = (20, "Score")
    peptide_intensity: Tuple[int, str] = (23, "Intensity (Peptide)")
    sequence: Tuple[int, str] = (8, "Sequence")
    confidence: Tuple[int, str] = (6, "Best Conf (Peptide)")
