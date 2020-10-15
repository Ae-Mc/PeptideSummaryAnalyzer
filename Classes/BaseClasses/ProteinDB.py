from typing import List, Dict
from os import listdir, path


class ProteinDB(dict):
    _necessaryColumns = ["Accession", "Unused"]

    def __init__(self, folder: str = None):
        if folder:
            self.LoadFromFolder(folder)

    def LoadFromFolder(self, folder: str) -> None:
        pass

    def LoadFromDict(self,
                     dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        pass

    def GetSortedTableNums(self):
        """Возвращает отсортированный список номеров таблиц

        Returns:
            Отсортированный список номеров таблиц, например:
            ["0.1", "0.2", ..., "10.1"]
        """
        return sorted(self, key=lambda x: float(x))

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
