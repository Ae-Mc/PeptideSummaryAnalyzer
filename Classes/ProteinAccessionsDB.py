from typing import Dict, List
from os import listdir, path
from decimal import Decimal
from Classes.Sequence import Sequence
from Classes.ReadTable import ReadTable
from Classes.Errors import ColumnNotFoundError
from Classes.ProteinAccession import ProteinAccession


class ProteinAccessionsDB(dict):
    """Хранит базу данных с Protein Accession в формате: {
        "номер таблицы": {
            "имя Accession": ProteinAccession
        }
    }
    """
    _necessaryColumns = ["Accession", "N", "Unused"]

    def __init__(self, folder: str = None):
        """См. LoadFromFolder"""
        if folder is not None:
            self.LoadFromFolder(folder)

    def LoadFromFolder(self, folder: str) -> None:
        """Загружает Protein таблицы из папки folder

        Args:
            folder: путь, в котором хранятся Protein таблицы
        """
        filenames = self._GetProteinFilenames(folder)
        dictionary: Dict[str, Dict[str, List[str]]] = {}
        for filename in filenames:
            tableNum = path.split(filename)[1].split('_')[0]
            dictionary[tableNum] = ReadTable(filename, unsafeFlag=True)
        self.LoadFromDict(dictionary)

    def LoadFromDict(self,
                     dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц

        Args:
            dictionary: словарь вида: {
                    "номер таблицы": {
                        "столбец": ["значение1", "значение2", ..., "значениеN"]
                    }
                }
        """
        for tableNum, table in dictionary.items():
            for columnName in self._necessaryColumns:
                if columnName not in table:
                    raise ColumnNotFoundError(columnName,
                                              f"{tableNum}_ProteinSummary.txt")
            curUnused: Decimal = Decimal(0)
            i = 0
            while i < len(table["Accession"]):
                if table["Accession"][i].startswith("RRRRR"):
                    break

                if Decimal(table["Unused"][i]) != Decimal(0):
                    curUnused = Decimal(table["Unused"][i])
                accession = ProteinAccession(table["Accession"][i],
                                             curUnused)
                if accession.name in self:
                    self[accession.name].unused = max(
                        self[accession.name].unused, accession.unused)
                    self[accession.name].occurences += 1
                else:
                    self[accession.name] = accession
                i += 1

    def GetRepresentative(self,
                          accessions: List[str],
                          seqDB: Dict[str, Sequence]) -> str:
        """Получает репрезентативный Accession для списка Accession

        Args:
            accessions: список имён Accession
            seqDB: база данных последовательностей Accession

        Returns:
            Имя репрезентативного Accession
        """
        representativeAccession = ProteinAccession("", Decimal(0), 0)
        for accession in sorted(accessions):
            if(
                self[accession].unused > representativeAccession.unused or
                (self[accession].unused == representativeAccession.unused and
                 self[accession].occurences >
                 representativeAccession.occurences) or
                (self[accession].unused == representativeAccession.unused and
                 self[accession].occurences ==
                 representativeAccession.occurences and
                 seqDB[accession].len >
                 seqDB[representativeAccession.name].len)):
                representativeAccession = self[accession]
        return representativeAccession.name

    @staticmethod
    def _GetProteinFilenames(folder: str) -> List[str]:
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
