from decimal import Decimal
from typing import List
from .BaseClasses.Table import Table
from .PeptideAccession import PeptideAccession
from .PeptideColumns import PeptideColumns


class PeptideTable(Table):
    """Считывает PeptideTable в список вида List[PeptideAccession]"""
    columns: PeptideColumns

    def __init__(self,
                 tableFilename: str = None,
                 unsafeFlag: bool = False,
                 columns: PeptideColumns = PeptideColumns()):
        self.columns = columns
        super().__init__(tableFilename, unsafeFlag)

    def Load(self, tableFilename) -> List[PeptideAccession]:
        super().Load(tableFilename)
        self.columns.TestColumnNames(self.pop(0))
        for i, line in enumerate(self):
            self[i] = PeptideAccession(
                name=line[self.columns.accession[0]],
                confidence=Decimal(line[self.columns.confidence[0]]),
                sc=Decimal(line[self.columns.sc[0]]),
                precursorSignal=(
                    Decimal(line[self.columns.precursorSignal[0]])
                    if len(line[self.columns.precursorSignal[0]].strip()) > 0
                    else Decimal(0)
                ),
                sequence=line[self.columns.sequence[0]]
            )
        return self

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return '[\n ' + "\n ".join([str(e) for e in self]) + '\n]'
