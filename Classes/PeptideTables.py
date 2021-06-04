from Classes.PeptideAccession import PeptideAccession
from Classes.Sequence import Sequence
from collections import defaultdict
from Classes.PeptideRow import PeptideRow
from Classes.RawPeptideTable import RawPeptideTable
from Classes.RawPeptideTables import RawPeptideTables
from typing import Dict, List, Tuple
from .PeptideTable import PeptideTable
from .ProteinPerTableList import ProteinPerTableList


class PeptideTables(dict):
    """Словарь вида {
        "номер таблицы": PeptideTable,
    }

    Attributes:
        rawPeptideTables: считанные построчно Peptide таблицы
        seqDB: база данных последовательностей Accession
    """

    rawPeptideTables: RawPeptideTables
    seqDB: Dict[str, Sequence]

    def __init__(
        self,
        rawPeptideTables: RawPeptideTables,
        seqDB: Dict[str, Sequence],
    ) -> None:
        """
        Args:
            rawPeptideTables: считанные построчно Peptide таблицы
            seqDB: база данных последовательностей Accession
        """

        self.rawPeptideTables = rawPeptideTables
        self.seqDB = seqDB
        self.ExtractFromRawPeptideTables()
        self.sortedTableNums = self.GetSortedTableNums()

    def ExtractFromRawPeptideTables(self) -> None:
        """Преобразовывает rawPeptideTables в PeptideTables"""
        tableNum: str
        table: RawPeptideTable
        countsPerTable = self.GetAccessionCountsPerTable()
        generalCounts = self.GetGeneralAccessionCounts(countsPerTable)
        for tableNum, table in self.rawPeptideTables.items():
            self[tableNum] = PeptideTable()
            row: PeptideRow
            print(tableNum)
            for row in table:
                representativeAccession = self.GetRepresentativeForRow(
                    row,
                    countsPerCurrentTable=countsPerTable[tableNum],
                    generalCounts=generalCounts,
                )
                self[tableNum].append(
                    PeptideAccession(
                        representativeAccession,
                        confidence=row.confidence,
                        sc=row.sc,
                        precursorSignal=row.precursorSignal,
                        sequence=row.sequence,
                    )
                )

    def GetAccessionCountsPerTable(self) -> Dict[str, Dict[str, int]]:
        countsPerTable = {}
        for tableNum, table in self.rawPeptideTables.items():
            countsPerTable[tableNum] = self.CountAccessionsInRawTable(table)
        return countsPerTable

    @staticmethod
    def CountAccessionsInRawTable(
        rawPeptideTable: RawPeptideTable,
    ) -> Dict[str, int]:
        counts: Dict[str, int] = defaultdict(lambda: 0)
        row: PeptideRow
        for row in rawPeptideTable:
            for accession in row.accessions:
                counts[accession] += 1
        return counts

    @staticmethod
    def GetGeneralAccessionCounts(
        countsPerTable: Dict[str, Dict[str, int]]
    ) -> Dict[str, int]:
        generalCounts: Dict[str, int] = defaultdict(lambda: 0)
        for table in countsPerTable.values():
            for accession, count in table.items():
                generalCounts[accession] += count
        return generalCounts

    def GetRepresentativeForRow(
        self,
        row: PeptideRow,
        countsPerCurrentTable: Dict[str, int],
        generalCounts: Dict[str, int],
    ) -> str:
        """Выбирает репрезентативный Accession для строки.

        Сначала учитывается количество появлений каждого Accession в текущей
        таблице, потом, если репрезентативный Accession не был выбран,
        учитывается учитывается количество появлений во всех таблицах,
        потом сравнивается длина последовательностей.

        Args:
            row: PeptideRow
            countsPerCurrentTable: количество появлений в текущей таблице
            generalCounts: количество появлений во всех таблицах

        Returns:
            str - репрезентативный Accession
        """
        representativeAccessions = self.GetRepresentativeForRowByCounts(
            row.accessions, countsPerCurrentTable
        )
        if len(representativeAccessions) == 1:
            return representativeAccessions[0]
        representativeAccessions = self.GetRepresentativeForRowByCounts(
            representativeAccessions, generalCounts
        )
        if len(representativeAccessions) == 1:
            return representativeAccessions[0]
        representativeAccessions = self.GetRepresentativeForRowBySequences(
            representativeAccessions
        )
        return representativeAccessions[0]

    @staticmethod
    def GetRepresentativeForRowByCounts(
        accessions: List[str], counts: Dict[str, int]
    ) -> List[str]:
        """Выбирает репрезентативный Accession на основе словаря с количеством
        появлений каждого Accession.

        Args:
            accessions: список Accession, для которого выбирается
                репрезентативный
            counts: словарь с количеством появлений каждого Accession

        Returns:
            None: если
        """
        representativeAccession: Tuple[str, int] = ("", 0)
        variants = []
        for accession in accessions:
            if counts[accession] > representativeAccession[1]:
                representativeAccession = (
                    accession,
                    counts[accession],
                )
                variants = [accession]
            elif counts[accession] == representativeAccession[1]:
                variants.append(accession)
        return variants

    def GetRepresentativeForRowBySequences(
        self, accessions: List[str]
    ) -> List[str]:
        representativeAccession: Tuple[str, int] = ("", 0)
        variants = []
        for accession in accessions:
            if self.seqDB[accession].len > representativeAccession[1]:
                representativeAccession = (
                    accession,
                    self.seqDB[accession].len,
                )
                variants = [accession]
            elif self.seqDB[accession].len == representativeAccession[1]:
                variants.append(accession)
        return variants

    def GetSortedTableNums(self) -> List[str]:
        """Получает отсортированный список номеров таблиц

        Returns:
            отсортированный список номеров таблиц
        """
        return sorted(self.keys(), key=lambda x: float(x))

    def ApplyProteinPerTableList(
        self, proteinPerTableList: ProteinPerTableList
    ) -> None:
        """Удаляет все Accession, отсутствующие в Protein таблицах

        Args:
            proteinPerTableList: список Protein Accession, распределённых по
                таблицам
        """
        for tableNum, table in self.items():
            i = 0
            while i < len(table):
                if table[i].name not in proteinPerTableList[tableNum]:
                    table.pop(i)
                    continue
                i += 1

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
