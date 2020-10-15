from dataclasses import dataclass
from .BaseClasses import BaseProteinAccession


@dataclass
class ProteinAccession(BaseProteinAccession):
    """Класс для хранения данных об Accession из файлов ProteinSummary

    Attributes:
        name: имя Accession
        unused: максимальный unused для данного Accession
        occurences: количество появлений Accession
    """
    occurences: int = 1

    def __str__(self):
        return (f"{self.name} (unused: {self.unused}, "
                f"occurences: {self.occurences})")

    def __repr__(self):
        return (f"{self.name} (unused: {self.unused}, "
                f"occurences: {self.occurences})")
