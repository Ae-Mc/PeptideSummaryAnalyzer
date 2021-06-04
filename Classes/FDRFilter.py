from Classes.Functions import IsReversed
from Classes.PeptideRow import PeptideRow
from Classes.RawPeptideTable import RawPeptideTable
from Classes.RawPeptideTables import RawPeptideTables


class FDRFilter:
    """Применяет FDR фильтр к RawPeptideTables

    Attributes:
        rawPeptideTables: считанные построчно Peptide таблицы
    """

    rawPeptideTables: RawPeptideTables

    def __init__(self, rawPeptideTables: RawPeptideTables) -> None:
        """
        Args:
            rawPeptideTables: считанные построчно Peptide таблицы
        """

        self.rawPeptideTables = rawPeptideTables

    def ApplyDefaultFilter(self):
        table: RawPeptideTable
        for table in self.rawPeptideTables.values():
            row: PeptideRow
            i = 0
            while i < len(table):
                row = table[i]
                self.RemoveReversedAccessionsFromRow(row)
                if len(row.accessions) == 0:
                    del table[i:]
                else:
                    i += 1

    @staticmethod
    def RemoveReversedAccessionsFromRow(row: PeptideRow) -> None:
        i = 0
        while i < len(row.accessions):
            if IsReversed(row.accessions[i]):
                row.accessions.pop(i)
            else:
                i += 1
