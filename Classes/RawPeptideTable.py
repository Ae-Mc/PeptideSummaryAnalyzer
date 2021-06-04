from Classes.PeptideRow import PeptideRow
from decimal import Decimal
from typing import List
from .BaseClasses.Table import Table
from .PeptideAccession import PeptideAccession
from .PeptideColumns import PeptideColumns


class RawPeptideTable(Table):
    """Считывает PeptideTable в список вида List[PeptideRow]."""

    columns: PeptideColumns

    def __init__(
        self,
        tableFilename: str = None,
        unsafeFlag: bool = False,
        columns: PeptideColumns = PeptideColumns(),
    ):
        """Initializes table settings and reads table from file"""
        self.columns = columns
        super().__init__(tableFilename, unsafeFlag)

    def Load(self, tableFilename) -> List[PeptideAccession]:
        super().Load(tableFilename)
        self.columns.TestColumnNames(self.pop(0))
        for i, line in enumerate(self):
            self[i] = PeptideRow(
                accessions=list(
                    map(
                        lambda x: x.strip(),
                        line[self.columns.accession[0]].split(";"),
                    )
                ),
                confidence=Decimal(line[self.columns.confidence[0]]),
                sc=Decimal(line[self.columns.sc[0]]),
                precursorSignal=(
                    Decimal(line[self.columns.precursorSignal[0]])
                    if len(line[self.columns.precursorSignal[0]].strip()) > 0
                    else Decimal(0)
                ),
                sequence=line[self.columns.sequence[0]],
            )
        return self

    def RemoveRowsWithAccessions(self, accessions: List[str]) -> None:
        i = 0
        while i < len(self):
            for accession in self[i].accessions:
                if accession in accessions:
                    self.pop(i)
                    break
            else:
                i += 1

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "[\n " + "\n ".join([str(e) for e in self]) + "\n]"
