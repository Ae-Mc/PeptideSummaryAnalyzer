from os import listdir
from typing import Dict, List
from .ReadTable import ReadTable
from .ColumnNames import ColumnNames
from .ProteinPerTableList import ProteinPerTableList


class PeptideTables(dict):
    """Словарь вида {
        "номер таблицы": {
            "Заголовок":  ["Значение1", "Значение2", ..., "ЗначениеN"],
        },
    }

    Attributes:
        columnNames: имена заголовков
    """
    columnNames: ColumnNames

    def __init__(self, columnNames: ColumnNames, inputDir: str = None) -> None:
        """
        Args:
            columnNames: имена заголовков
            inputDir: путь, из которого считываются таблицы
        """

        self.SetColumnNames(columnNames)

        if inputDir is not None:
            self.ReadPeptideSummaries(inputDir)
            self.sortedTableNums = self.GetSortedTableNums()
            self.RemoveReversedAccessions()
            self.RemoveExcessAccessions()

    def ReadPeptideSummaries(self, inputDir: str) -> None:
        """Считывает все PeptideSummary файлы в словарь

        Args:
            inputDir: путь, из которого считываются таблицы
        """
        for filename in listdir(inputDir):
            if "Peptide" in filename:
                tableNum = filename.split('_')[0]
                self[tableNum] = (
                    ReadTable(inputDir + '/' + filename))
                curTableColumnNames = list(self[tableNum].keys())
                for columnName in curTableColumnNames:
                    if columnName not in self.columnNames.GetColumnNamesList():
                        del self[tableNum][columnName]

    def GetSortedTableNums(self) -> List[str]:
        """Получает отсортированный список номеров таблиц

        Returns:
            отсортированный список номеров таблиц
        """
        return sorted(self.keys(), key=lambda x: float(x))

    def RemoveReversedAccessions(self) -> None:
        """Удаляет перевёрнутые Accession

        Перевёрнутые Accession - это Accession, начинающиеся с RRRRR
        """
        for peptideTable in self.values():
            i = 0
            while i < len(peptideTable[self.columnNames.accession]):
                if peptideTable[
                        self.columnNames.accession][i].startswith("RRRRR"):
                    break
                i += 1

            if i < len(peptideTable[self.columnNames.accession]):
                for column in peptideTable.values():
                    del column[i:]

    def RemoveExcessAccessions(self) -> None:
        """Удаляет все лишние имена, которые идут после ; в имени
        Accession
        """
        for table in self.values():
            table[self.columnNames.accession] = [
                accession.split(';')[0] for accession in table[
                    self.columnNames.accession]]

    def ApplyProteinPerTableList(
            self, proteinPerTableList: ProteinPerTableList) -> None:
        """Удаляет все Accession, отсутствующие в Protein таблицах

        Args:
            proteinPerTableList: список Protein Accession, распределённых по
                таблицам
        """
        for tableNum, table in self.items():
            i = 0
            while i < len(table[self.columnNames.accession]):
                if(table[self.columnNames.accession][i] not in
                   proteinPerTableList[tableNum]):
                    self.RemoveRow(tableNum, i)
                    i -= 1
                i += 1

    def ApplyProteinReplacements(
            self, proteinReplacements: Dict[str, Dict[str, str]]) -> None:
        """Применяет замены, полученные из Protein таблиц

        Args:
            proteinReplacements: словарь замен вида: {
                    "номер таблицы": {
                        "что заменяем": "на что заменяем"
                    }
                }
        """
        for tableNum, table in self.items():
            tableReplacements = proteinReplacements[tableNum]
            for i in range(0, len(table[self.columnNames.accession])):
                if table[self.columnNames.accession][i] in tableReplacements:
                    table[self.columnNames.accession][i] = (tableReplacements[
                        table[self.columnNames.accession][i]])

    def SetColumnNames(self, columnNames: ColumnNames):
        self.columnNames = columnNames

    def RemoveRow(self, tableNum: str, rowNum: int) -> None:
        columns = [column for column in self[tableNum]]
        for column in columns:
            del self[tableNum][column][rowNum]
