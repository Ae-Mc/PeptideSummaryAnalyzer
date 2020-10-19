from typing import Tuple
from .PeptideColumns import PeptideColumns, dataclass


@dataclass
class PeptideColumns5(PeptideColumns):

    """См. PeptideColumns"""

    sc: Tuple[int, str] = (22, "Sc")
    precursorSignal: Tuple[int, str] = (25, "Intensity (Peptide)")
