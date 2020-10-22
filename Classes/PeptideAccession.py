from dataclasses import dataclass
from decimal import Decimal


@dataclass
class PeptideAccession:
    name: str = ""
    unused: Decimal = Decimal(0)
    confidence: Decimal = Decimal(0)
    contribution: Decimal = Decimal(0)
    sc: Decimal = Decimal(0)
    precursorSignal: Decimal = Decimal(0)
    sequence: str = ""
