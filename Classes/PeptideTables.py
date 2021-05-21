from os import listdir
from typing import Dict, List
from .PeptideColumns import PeptideColumns
from .PeptideTable import PeptideTable
from .ProteinPerTableList import ProteinPerTableList


class PeptideTables(dict):
    """Словарь вида {
        "номер таблицы": PeptideTable,
    }

    Attributes:
        columnNames: имена заголовков
    """
    columnNames: PeptideColumns

    def __init__(self,
                 columnNames: PeptideColumns,
                 inputDir: str = None) -> None:
        """
        Args:
            columnNames: имена заголовков
            inputDir: путь, из которого считываются таблицы
        """

        self.columnNames = columnNames

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
                self[tableNum] = PeptideTable(inputDir + '/' + filename,
                                              unsafeFlag=True,
                                              columns=self.columnNames)

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
        peptideTable: PeptideTable
        for peptideTable in self.values():
            i = 0
            while i < len(peptideTable):
                if peptideTable[i].name.startswith("RRRRR"):
                    break
                i += 1

            del peptideTable[i:]

    def RemoveExcessAccessions(self) -> None:
        """Удаляет все лишние имена, которые идут после ; в имени
        Accession
        """
        for table in self.values():
            for i in range(len(table)):
                table[i].name = table[i].name.split(';')[0]

    def ApplyProteinPerTableList(
            self, proteinPerTableList: ProteinPerTableList) -> None:
        """Удаляет все Accession, отсутствующие в Protein таблицах

        Args:
            proteinPerTableList: список Protein Accession, распределённых по
                таблицам
        """
        for tableNum, table in self.items():
            i = 0
            while i < len(table):
                if(table[i].name not in proteinPerTableList[tableNum]):
                    table.pop(i)
                    continue
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
            for i in range(len(table)):
                if table[i].name in tableReplacements:
                    table[i].name = tableReplacements[table[i].name]
