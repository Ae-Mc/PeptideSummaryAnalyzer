from decimal import Decimal
from typing import List
from .BaseClasses.Table import Table
from .BaseClasses.BaseProteinAccession import BaseProteinAccession
from .ProteinColumns import ProteinColumns


class ProteinTable(Table):
    columns: ProteinColumns

    def __init__(self,
                 tableFilename: str = None,
                 unsafeFlag: bool = False,
                 columns: ProteinColumns = ProteinColumns()):
        self.columns = columns
        super().__init__(tableFilename, unsafeFlag)

    def Load(self, tableFilename: str) -> List[BaseProteinAccession]:
        super().Load(tableFilename)
        self.columns.TestColumnNames(self.pop(0))
        for i, line in enumerate(self):
            self[i] = BaseProteinAccession(line[self.columns.accession[0]],
                                           Decimal(
                                               line[self.columns.unused[0]]))
        return self
