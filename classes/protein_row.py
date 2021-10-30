"""См. класс ProteinRow."""

from dataclasses import dataclass
from decimal import Decimal


@dataclass
class ProteinRow:
    """Хранит информацию о строках из ProteinSummary файлов

    Attributes:
        n (int): значение N строки
        accession (str): Accession
        unused (Decimal): значение unused для строки
    """

    n: int  # pylint: disable=invalid-name
    accession: str
    unused: Decimal
