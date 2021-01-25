from decimal import Decimal
from typing import List
from .BaseClasses.Table import Table
from .PeptideAccession import PeptideAccession
from .PeptideColumns import PeptideColumns


class PeptideTable(Table):
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
                unused=Decimal(line[self.columns.unused[0]]),
                confidence=Decimal(line[self.columns.confidence[0]]),
                sc=Decimal(line[self.columns.sc[0]]),
                precursorSignal=Decimal(line[self.columns.precursorSignal[0]]),
                sequence=line[self.columns.sequence[0]]
            )
        return self
