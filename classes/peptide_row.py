"""См. класс PeptideRow."""

from dataclasses import dataclass
from decimal import Decimal
from typing import List


@dataclass
class PeptideRow:
    """Хранит информацию о группе из DistinctPeptideSummary файлов

    Attributes:
        accessions (List[str]): Список имён Accession
        confidence (Decimal): значение confidence для группы
        score (Decimal): значение sc для группы
        peptide_intensity (Decimal): значение Peptide Intensity для группы
        sequence (str): строка последовательности группы
    """

    accessions: List[str]
    confidence: Decimal
    score: Decimal
    peptide_intensity: Decimal
    sequence: str
