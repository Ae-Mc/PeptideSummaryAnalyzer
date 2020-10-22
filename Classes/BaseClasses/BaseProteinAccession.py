from decimal import Decimal
from dataclasses import dataclass


@dataclass
class BaseProteinAccession:
    """Класс для хранения данных об Accession из файлов ProteinSummary

    Attributes:
        name: имя Accession
        unused: максимальный unused для данного Accession
    """
    name: str
    unused: Decimal

    def __str__(self):
        return f"{self.name} (unused: {self.unused})"

    def __repr__(self):
        return f"{self.name} (unused: {self.unused})"
