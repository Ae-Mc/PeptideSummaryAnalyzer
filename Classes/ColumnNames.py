from dataclasses import dataclass


@dataclass
class ColumnNames:

    """Хранит названия столбцов для таблиц Peptide

    Attributes:
        accession: имя столбца accession
        precursorSignal: имя столбца precursorSignal
        sc: имя столбца sc
        unused: имя столбца unused
        sequence: имя столбца sequence
        confidence: имя столбца confidence
        contribution: имя столбца contribution
    """

    accession: str = "Accessions"
    precursorSignal: str = "PrecursorSignal"
    sc: str = "Sc"
    unused: str = "Unused"
    sequence: str = "Sequence"
    confidence: str = "Conf"
    contribution: str = "Contrib"

    def GetColumnNamesList(self):
        return [self.accession,
                self.precursorSignal,
                self.sc,
                self.unused,
                self.sequence,
                self.confidence,
                self.contribution]
