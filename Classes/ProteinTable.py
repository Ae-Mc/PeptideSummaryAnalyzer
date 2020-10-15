from .BaseClasses import Table, BaseProteinAccession
from .ProteinColumns import ProteinColumns
from typing import List


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
        for i, line in enumerate(self):
            self[i] = BaseProteinAccession(
                line[self.columns.accession], line[self.columns.unused])
        return self
