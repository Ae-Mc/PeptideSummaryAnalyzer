from os import listdir
from typing import Dict, List
from Classes.ReadTable import ReadTable
from Classes.ProteinTables import ProteinTables
from Classes.ColumnNames import ColumnNames

""" peptideTables - словарь вида
    {
        "1.1": {
            "Заголовок 1":  ["Значение 1", "Значение 2", ..., "Значение n1"],
            "Заголовок 2":  ["Значение 1", "Значение 2", ..., "Значение n1"],
            ..............................................................
            "Заголовок N1": ["Значение 1", "Значение 2", ..., "Значение n1"]
        },
        "1.2": {
            "Заголовок 1":  ["Значение 1", "Значение 2", ..., "Значение n2"],
            "Заголовок 2":  ["Значение 1", "Значение 2", ..., "Значение n2"],
            ..............................................................
            "Заголовок N2": ["Значение 1", "Значение 2", ..., "Значение n2"]
        },
        ...,
        "I-ый файл": {
            "Заголовок 1":  ["Значение 1", "Значение 2", ..., "Значение nI"],
            "Заголовок 2":  ["Значение 1", "Значение 2", ..., "Значение nI"],
            ..............................................................
            "Заголовок NI": ["Значение 1", "Значение 2", ..., "Значение nI"]
        }
    }
"""


class PeptideTables:

    peptideTables: Dict[str, Dict[str, List[str]]]
    columnNames: ColumnNames

    def __init__(self,
                 inputDir: str = None,
                 columnNames: ColumnNames = None):

        if columnNames is None:
            self.SetColumnNames(ColumnNames())
        else:
            self.SetColumnNames(columnNames)

        if inputDir is not None:
            self.ReadPeptideSummaries(inputDir)
            self.sortedTableNums = self.GetSortedTableNums()
            self.RemoveReversedAccessions()
            self.RemoveExcessAccessions()

    def ReadPeptideSummaries(self, inputDir: str) -> None:
        """ Считывание всех PeptideSummary файлов в словарь """

        self.peptideTables = {}
        for filename in listdir(inputDir):
            if "Peptide" in filename:
                tableNum = filename.split('_')[0]
                self.peptideTables[tableNum] = (
                    ReadTable(inputDir + '/' + filename))

    def GetSortedTableNums(self) -> List[str]:
        return sorted(self.peptideTables.keys(), key=lambda x: float(x))

    def RemoveReversedAccessions(self):

        for peptideTable in self.peptideTables.values():
            i = 0
            while i < len(peptideTable[self.columnNames.accession]):
                if peptideTable[
                        self.columnNames.accession][i].startswith("RRRRR"):
                    break
                i += 1

            if i < len(peptideTable[self.columnNames.accession]):
                for column in peptideTable.values():
                    del column[i:]

    def RemoveExcessAccessions(self):

        for table in self.peptideTables.values():
            table[self.columnNames.accession] = [
                accession.split(';')[0] for accession in table[
                    self.columnNames.accession]]

    def ApplyProteinReplacements(self, proteinTables: ProteinTables):

        for tableName, table in self.peptideTables.items():
            tableReplacements = proteinTables.proteinReplacements[tableName]
            for i in range(0, len(table[self.columnNames.accession])):
                if table[self.columnNames.accession][i] in tableReplacements:
                    table[self.columnNames.accession][i] = (tableReplacements[
                        table[self.columnNames.accession][i]])

    def SetColumnNames(self, columnNames: ColumnNames):
        self.columnNames = columnNames
