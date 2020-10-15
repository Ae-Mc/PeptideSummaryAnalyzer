from typing import Dict, List
from os import path
from decimal import Decimal
from Classes.ProteinAccessionsDB import ProteinAccessionsDB
from Classes.ProteinGroup import ProteinGroup
from Classes.ReadTable import ReadTable
from Classes.Errors import ColumnNotFoundError, AccessionNotFoundError
from Classes.Sequence import Sequence


class ProteinGroupsDB(dict):
    necessaryColumns = ["Accession", "Unused"]
    sortedTableNums: List[str]
    proteinAccessionsDB: ProteinAccessionsDB
    seqDB: Dict[str, Sequence]
    skipReversedIfSecondary: bool

    def __init__(self,
                 seqDB: Dict[str, Sequence],
                 skipReversedIfSecondary: bool = False,
                 folder: str = None) -> None:
        """См. LoadFromFolder

        Args:
            proteinAccessionsDB: база данных с Protein таблицами
            seqDB: база данных последовательностей Accession
            folder: путь, в котором хранятся Protein таблицы
        """
        self.skipReversedIfSecondary = skipReversedIfSecondary
        self.seqDB = seqDB
        if folder:
            self.LoadFromFolder(folder)

    def LoadFromFolder(self, folder: str) -> None:
        """Загружает Protein таблицы из папки folder

        Args:
            folder: путь, в котором хранятся Protein таблицы
        """
        self.proteinAccessionsDB = ProteinAccessionsDB(
            self.skipReversedIfSecondary, folder)
        filenames = ProteinAccessionsDB._GetProteinFilenames(folder)
        dictionary: Dict[str, Dict[str, List[str]]] = {}
        for filename in filenames:
            tableNum = path.split(filename)[1].split('_')[0]
            dictionary[tableNum] = ReadTable(filename, unsafeFlag=True)
        self.LoadFromDict(dictionary)

    def LoadFromDict(self,
                     dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц и создаёт из них группы, а также вычисляет для
        них репрезентативные Accession

        Args:
            dictionary: словарь вида: {
                    "номер таблицы": {
                        "столбец": ["значение1", "значение2", ..., "значениеN"]
                    }
                }
        """
        self.GetAccessionsGroups(dictionary)
        self.CalculateRepresentatives()

    def GetAccessionsGroups(
            self, dictionary: Dict[str, Dict[str, List[str]]]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц и создаёт из них группы

        Args:
            dictionary: словарь вида: {
                    "номер таблицы": {
                        "столбец": ["значение1", "значение2", ..., "значениеN"]
                    }
                }
        """
        self.sortedTableNums = sorted([*dictionary], key=lambda x: float(x))
        for tableNum, table in dictionary.items():
            self[tableNum] = []
            self.TestNeccessaryColumnNames(table, tableNum)
            self.GetAccessionsGroupsFromTable(table, tableNum)

    def TestNeccessaryColumnNames(
            self, table: Dict[str, List[str]], tableNum: str) -> None:
        for columnName in self.necessaryColumns:
            if columnName not in table:
                raise ColumnNotFoundError(columnName,
                                          f"{tableNum}_ProteinSummary.txt")

    def GetAccessionsGroupsFromTable(
            self, table: Dict[str, List[str]], tableNum: str) -> None:
        if len(table["Accession"]):
            curGroup = ProteinGroup(
                Decimal(table["Unused"][0]), [table["Accession"][0]])
            reversedFound = False
            i = 1
            while i < len(table["Accession"]):
                if (table["Accession"][i].startswith("RRRRR")):
                    if self.TestIsSecondaryReversed(table["Accession"], i):
                        return
                    if reversedFound:
                        return
                    reversedFound = True
                    i += 1
                    continue

                if Decimal(table["Unused"][i]) != Decimal(0):
                    curGroup.accessions = sorted(curGroup.accessions)
                    self[tableNum].append(curGroup)
                    curGroup = ProteinGroup(
                        Decimal(table["Unused"][i]), [])
                curGroup.accessions.append(table["Accession"][i])
                i += 1
            curGroup.accessions = sorted(curGroup.accessions)
            self[tableNum].append(curGroup)

    def TestIsSecondaryReversed(
            self, accessionColumn: List[str], i: int) -> bool:
        return not (not self.skipReversedIfSecondary
                    or i + 1 == len(accessionColumn)
                    or accessionColumn[i + 1].startswith("RRRRR"))

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
