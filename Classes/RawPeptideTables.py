from Classes.RawPeptideTable import RawPeptideTable
from os import listdir
from typing import Dict, List
from .PeptideColumns import PeptideColumns


class RawPeptideTables(dict):
    """Словарь вида {
        "номер таблицы": RawPeptideTable,
    }

    Attributes:
        columnNames: имена заголовков
    """

    columnNames: PeptideColumns

    def __init__(
        self, columnNames: PeptideColumns, inputDir: str = None
    ) -> None:
        """
        Args:
            columnNames: имена заголовков
            inputDir: путь, из которого считываются таблицы
        """

        self.columnNames = columnNames

        if inputDir is not None:
            self.ReadPeptideSummaries(inputDir)
            self.sortedTableNums = self.GetSortedTableNums()

    def ReadPeptideSummaries(self, inputDir: str) -> None:
        """Считывает все PeptideSummary файлы в словарь

        Args:
            inputDir: путь, из которого считываются таблицы
        """
        for filename in listdir(inputDir):
            if "Peptide" in filename:
                tableNum = filename.split("_")[0]
                self[tableNum] = RawPeptideTable(
                    inputDir + "/" + filename,
                    unsafeFlag=True,
                    columns=self.columnNames,
                )

    def GetSortedTableNums(self) -> List[str]:
        """Получает отсортированный список номеров таблиц

        Returns:
            отсортированный список номеров таблиц
        """
        return sorted(self.keys(), key=lambda x: float(x))

    def ApplyProteinReplacements(
        self, proteinReplacements: Dict[str, Dict[str, str]]
    ) -> None:
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
