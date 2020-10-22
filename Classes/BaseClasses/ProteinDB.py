from abc import ABC, abstractmethod
from typing import Dict, List
from os import path, listdir
from ..ProteinTable import ProteinTable


class ProteinDB(dict, ABC):
    _necessaryColumns = ["Accession", "Unused"]

    def __init__(self, proteinTables: Dict[str, ProteinTable] = None):
        if proteinTables:
            self.LoadFromTables(proteinTables)

    def LoadFromFolder(self, folder: str) -> None:
        """Загружает Protein таблицы из папки folder

        Args:
            folder: путь, в котором хранятся Protein таблицы
        """
        self.LoadFromTables(self.GetProteinTables(folder))

    @abstractmethod
    def LoadFromTables(self, dictionary: Dict[str, ProteinTable]) -> None:
        pass

    def GetSortedTableNums(self):
        """Возвращает отсортированный список номеров таблиц

        Returns:
            Отсортированный список номеров таблиц, например:
            ["0.1", "0.2", ..., "10.1"]
        """
        return sorted(self, key=lambda x: float(x))

    @staticmethod
    def GetProteinTables(
            folder: str,
            skipReversedIfSecondary: bool = False) -> Dict[str, ProteinTable]:
        """Загружает все Protein файлы по пути folder в классы ProteinTable

        Args:
            folder: путь к папки с Protein файлами

        Returns:
            Словарь вида: {
                "Номер таблицы": ProteinTable
            }
        """
        filenames = ProteinDB.GetProteinFilenames(folder)
        tables: Dict[str, ProteinTable] = {}
        for filename in filenames:
            tableNum = path.split(filename)[1].split('_')[0]
            tables[tableNum] = ProteinTable(
                tableFilename=filename,
                skipReversedIfSecondary=skipReversedIfSecondary,
                unsafeFlag=True)
        return tables

    @staticmethod
    def GetProteinFilenames(folder: str) -> List[str]:
        """Получает имена всех Protein файлов по пути folder

        Args:
            folder: путь к папки с Protein файлами

        Returns:
            Список имён Protein файлов по пути folder
        """
        filenames: List[str] = []
        for filename in listdir(folder):
            if filename.endswith("ProteinSummary.txt"):
                filenames.append(path.join(folder, filename))
        return filenames
