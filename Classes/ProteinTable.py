from decimal import Decimal
from typing import List
from .BaseClasses.Table import Table
from .BaseClasses.BaseProteinAccession import BaseProteinAccession
from .ProteinColumns import ProteinColumns


class ProteinTable(Table):
    columns: ProteinColumns
    skipReversedIfSecondary: bool

    def __init__(self,
                 tableFilename: str = None,
                 skipReversedIfSecondary: bool = False,
                 unsafeFlag: bool = False,
                 columns: ProteinColumns = ProteinColumns()):
        self.columns = columns
        self.skipReversedIfSecondary = skipReversedIfSecondary
        super().__init__(tableFilename, unsafeFlag)

    def Load(self, tableFilename: str) -> List[BaseProteinAccession]:
        super().Load(tableFilename)
        self.columns.TestColumnNames(self.pop(0))
        if self.skipReversedIfSecondary:
            reversedFound = False
            i = 0
            while i < len(self):
                self[i] = BaseProteinAccession(
                    self[i][self.columns.accession[0]], Decimal(
                        self[i][self.columns.unused[0]]))
                if(self[i].name.startswith("RRRRR")):
                    if reversedFound:
                        break
                    if(i + 1 == len(self)
                       or self[i + 1][self.columns.accession[0]].startswith(
                            "RRRRR")):
                        break
                    reversedFound = True
                    self.pop(i)
                    continue
                i += 1
            else:
                return self
            del self[i:]
        else:
            for i, line in enumerate(self):
                self[i] = BaseProteinAccession(
                    line[self.columns.accession[0]], Decimal(
                        line[self.columns.unused[0]]))
                if(self[i].name.startswith("RRRRR")):
                    break
            else:
                return self
            del self[i:]
        return self
