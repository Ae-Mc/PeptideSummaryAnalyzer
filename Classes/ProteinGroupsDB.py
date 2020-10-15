from typing import Dict, List
from os import path
from decimal import Decimal
from .ProteinAccessionsDB import ProteinAccessionsDB
from .ProteinGroup import ProteinGroup
from .ReadTable import ReadTable
from .Errors import ColumnNotFoundError, AccessionNotFoundError
from .Sequence import Sequence
from .BaseClasses import ProteinDB


class ProteinGroupsDB(ProteinDB):
    """Хранит Protein группы, разбитые по таблицам

    Attributes:
        proteinAccessionsDB: база данных accession, не разбитая по группам
        seqDB: база данных последовательностей Accession
    """
    proteinAccessionsDB: ProteinAccessionsDB
    seqDB: Dict[str, Sequence]

    def __init__(self,
                 proteinAccessionsDB: ProteinAccessionsDB,
                 seqDB: Dict[str, Sequence],
                 folder: str) -> None:
        """См. LoadFromFolder

        Args:
            proteinAccessionsDB: база данных с Protein таблицами
            seqDB: база данных последовательностей Accession
            folder: путь, в котором хранятся Protein таблицы
        """
        self.proteinAccessionsDB = proteinAccessionsDB
        self.seqDB = seqDB
        self.LoadFromFolder(folder)

    def LoadFromFolder(self, folder: str) -> None:
        """Загружает Protein таблицы из папки folder

        Args:
            folder: путь, в котором хранятся Protein таблицы
        """
        filenames = ProteinDB.GetProteinFilenames(folder)
        dictionary: Dict[str, Dict[str, List[str]]] = {}
        for filename in filenames:
            tableNum = path.split(filename)[1].split('_')[0]
            dictionary[tableNum] = ReadTable(filename, unsafeFlag=True)
        self.LoadFromDict(dictionary)

    def LoadFromDict(self,
                     dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц и создаёт из них группы

        Args:
            dictionary: словарь вида: {
                    "номер таблицы": {
                        "столбец": ["значение1", "значение2", ..., "значениеN"]
                    }
                }
        """
        for tableNum, table in dictionary.items():
            self[tableNum] = []
            for columnName in self._necessaryColumns:
                if columnName not in table:
                    raise ColumnNotFoundError(
                        columnName, f"{tableNum}_ProteinSummary.txt")
            if len(table["Accession"]):
                curGroup = ProteinGroup(
                    Decimal(table["Unused"][0]), [table["Accession"][0]])
                i = 1
                while i < len(table["Accession"]):
                    if table["Accession"][i].startswith("RRRRR"):
                        break

                    if Decimal(table["Unused"][i]) != Decimal(0):
                        curGroup.accessions = sorted(curGroup.accessions)
                        self[tableNum].append(curGroup)
                        curGroup = ProteinGroup(
                            Decimal(table["Unused"][i]), [])
                    curGroup.accessions.append(table["Accession"][i])
                    i += 1
                self[tableNum].append(curGroup)
        self.CalculateRepresentatives()

    def CalculateRepresentatives(self):
        """Подсчитывает репрезентативные Accession для всех групп"""
        for tableNum, table in self.items():
            for group in table:
                group.representativeAccession = (
                    self.proteinAccessionsDB.GetRepresentative(
                        group.accessions, self.seqDB))

    def GetReplacementsPerTable(self) -> Dict[str, Dict[str, str]]:
        """Создаёт словарь замен

        Returns:
            словарь замен вида: {
                "номер таблицы": {
                    "что заменять": "на что заменять"
                }
            }
        """
        replacements: Dict[str, Dict[str, str]] = {}
        for tableNum, groups in self.items():
            replacements[tableNum] = {}
            group: ProteinGroup
            for group in groups:
                if group.representativeAccession is None:
                    raise AccessionNotFoundError(
                        "Representative accession for group"
                        f"{group.accessions} not found")
                for accession in group.accessions:
                    replacements[tableNum][accession] = (
                        group.representativeAccession)
        return replacements
