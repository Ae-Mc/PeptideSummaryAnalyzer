from os import listdir
from typing import Dict, List
from .ReadTable import ReadTable
from .ColumnNames import ColumnNames
from .ProteinPerTableList import ProteinPerTableList


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
                 columnNames: ColumnNames,
                 inputDir: str = None) -> None:

        self.SetColumnNames(columnNames)

        if inputDir is not None:
            self.ReadPeptideSummaries(inputDir)
            self.sortedTableNums = self.GetSortedTableNums()
            self.RemoveReversedAccessions()
            self.RemoveExcessAccessions()

    def ReadPeptideSummaries(self, inputDir: str) -> None:
        """ Считывание всех PeptideSummary файлов в словарь """

        self.peptideTables: Dict[str, Dict[str, List[str]]] = {}
        for filename in listdir(inputDir):
            if "Peptide" in filename:
                tableNum = filename.split('_')[0]
                self.peptideTables[tableNum] = (
                    ReadTable(inputDir + '/' + filename))
                curTableColumnNames = list(self.peptideTables[tableNum].keys())
                for columnName in curTableColumnNames:
                    if columnName not in self.columnNames.GetColumnNamesList():
                        del self.peptideTables[tableNum][columnName]

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

    def ApplyProteinPerTableList(
            self, proteinPerTableList: ProteinPerTableList) -> None:
        for tableNum, table in self.peptideTables.items():
            i = 0
            while i < len(table[self.columnNames.accession]):
                if(table[self.columnNames.accession][i] not in
                   proteinPerTableList[tableNum]):
                    self.RemoveRow(tableNum, i)
                    i -= 1
                i += 1

    def ApplyProteinReplacements(
            self, proteinReplacements: Dict[str, Dict[str, str]]):

        for tableNum, table in self.peptideTables.items():
            tableReplacements = proteinReplacements[tableNum]
            for i in range(0, len(table[self.columnNames.accession])):
                if table[self.columnNames.accession][i] in tableReplacements:
                    table[self.columnNames.accession][i] = (tableReplacements[
                        table[self.columnNames.accession][i]])

    def SetColumnNames(self, columnNames: ColumnNames):
        self.columnNames = columnNames

    def RemoveRow(self, tableNum: str, rowNum: int) -> None:
        columns = [column for column in self.peptideTables[tableNum]]
        for column in columns:
            del self.peptideTables[tableNum][column][rowNum]
