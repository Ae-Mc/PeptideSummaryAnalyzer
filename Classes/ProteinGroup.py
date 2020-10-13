from typing import List, Optional
from decimal import Decimal
from dataclasses import dataclass, field


@dataclass
class ProteinGroup:
    """Хранит Protein группу

    Attributes:
        unused: значение unused для всей группы
        accessions: список Accession
        representativeAccession: репрезентативный Accession
    """
    unused: Decimal = Decimal(0)
    accessions: List[str] = field(default_factory=list)
    representativeAccession: Optional[str] = None

    def __repr__(self):
        return (f"{self.representativeAccession} "
                f"({self.unused}): {self.accessions}")

    def __str__(self):
        return (f"{self.representativeAccession} "
                f"({self.unused}): {self.accessions}")
