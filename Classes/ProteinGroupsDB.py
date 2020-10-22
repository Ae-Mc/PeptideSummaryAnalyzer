from typing import Dict
from decimal import Decimal
from .ProteinAccessionsDB import ProteinAccessionsDB
from .ProteinGroup import ProteinGroup
from .ProteinTable import ProteinTable
from .Errors import AccessionNotFoundError
from .Sequence import Sequence
from .BaseClasses.ProteinDB import ProteinDB


class ProteinGroupsDB(ProteinDB):
    """Хранит Protein группы, разбитые по таблицам

    Attributes:
        proteinAccessionsDB: база данных accession, не разбитая по группам
        seqDB: база данных последовательностей Accession
    """
    proteinAccessionsDB: ProteinAccessionsDB
    seqDB: Dict[str, Sequence]

    def __init__(self,
                 seqDB: Dict[str, Sequence],
                 proteinTables: Dict[str, ProteinTable] = None) -> None:
        """См. LoadFromFolder

        Args:
            proteinAccessionsDB: база данных с Protein таблицами
            seqDB: база данных последовательностей Accession
            proteinTables: словарь ProteinTable вида: {
                    "Номер таблицы": ProteinTable
                }
        """
        self.seqDB = seqDB
        if proteinTables:
            self.LoadFromTables(proteinTables)

    def LoadFromTables(self, proteinTables: Dict[str, ProteinTable]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц и создаёт из них группы, а также вычисляет для
        них репрезентативные Accession

        Args:
            proteinTables: словарь вида: {
                    "номер таблицы": ProteinTable
                }
        """
        self.GetAccessionsGroups(proteinTables)
        self.CalculateRepresentatives(proteinTables)

    def GetAccessionsGroups(
            self, proteinTables: Dict[str, ProteinTable]) -> None:
        """Загружает Protein таблицы из словаря, полученного в результате
        чтения Protein таблиц и создаёт из них группы

        Args:
            proteinTables: словарь вида: {
                    "номер таблицы": ProteinTable
                }
        """
        for tableNum, table in proteinTables.items():
            self[tableNum] = []
            if len(table):
                curGroup = ProteinGroup(table[0].unused, [table[0].name])
                i = 1
                while i < len(table):
                    if table[i].unused != Decimal(0):
                        curGroup.accessions = sorted(curGroup.accessions)
                        self[tableNum].append(curGroup)
                        curGroup = ProteinGroup(table[i].unused, [])
                    curGroup.accessions.append(table[i].name)
                    i += 1
                curGroup.accessions = sorted(curGroup.accessions)
                self[tableNum].append(curGroup)

    def CalculateRepresentatives(self, proteinTables: Dict[str, ProteinTable]):
        """Подсчитывает репрезентативные Accession для всех групп

        Args:
            proteinTables: словарь вида: {
                    "Номер таблицы": ProteinTable
                }
        """
        proteinAccessionsDB = ProteinAccessionsDB(proteinTables)
        for tableNum, table in self.items():
            for group in table:
                group.representativeAccession = (
                    proteinAccessionsDB.GetRepresentative(group.accessions,
                                                          self.seqDB))

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
