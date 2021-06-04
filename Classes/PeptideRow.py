from decimal import Decimal
from typing import List


class PeptideRow:
    """Хранит информацию о группе из DistinctPeptideSummary файлов

    Attributes:
        accessions (List[str]): Список имён Accession
        confidence (Decimal): значение confidence для группы
        sc (Decimal): значение sc для группы
        precursorSignal (Decimal): значение precursorSignal для группы
        sequence (str): строка последовательности группы
    """

    accessions: List[str]
    confidence: Decimal
    sc: Decimal
    precursorSignal: Decimal
    sequence: str

    def __init__(
        self,
        accessions: List[str],
        confidence: Decimal,
        sc: Decimal,
        precursorSignal: Decimal,
        sequence: str,
    ):
        """Sets values for class instance."""

        self.accessions = accessions
        self.confidence = confidence
        self.sc = sc
        self.precursorSignal = precursorSignal
        self.sequence = sequence
