from typing import Dict, List
from Classes.Accession import Accession
from Classes.Sequence import Sequence
from Classes.PeptideTables import PeptideTables
from Classes.ColumnNames import ColumnNames
from decimal import Decimal


class AccessionTables:

    accessionsPerTable: Dict[str, Dict[str, Accession]]
    sortedTableNums: List[str]
    columnNames: ColumnNames

    def __init__(self,
                 seqDB: Dict[str, Sequence],
                 peptideTables: PeptideTables,
                 columnNames: ColumnNames = None) -> None:
        if columnNames is None:
            self.columnNames = ColumnNames()
        else:
            self.columnNames = columnNames
        self.proteinReplacements = None
        self.GetAccessionsPerTable(seqDB, peptideTables)

    def GetAccessionsPerTable(
            self,
            seqences: Dict[str, Sequence],
            peptideTables: PeptideTables) -> None:
        """ Получаем суммы значений Sc, Precursor Signal и сумму длинн
        последовательностей для каждого Accession, а также нормализованные
        значения Precursor Signal и Sc для каждого файла"""

        self.accessionsPerTable: Dict[str, Dict[str, Accession]] = {}
        for tableNum in peptideTables.peptideTables:
            self.accessionsPerTable[tableNum] = (
                self.GetAccessionsFromTable(
                    peptideTables.peptideTables[tableNum]))
            self.CalculateNormParamsForAccessions(
                self.accessionsPerTable[tableNum], seqences)

    def GetAccessionsFromTable(
            self,
            peptideTable: Dict[str, List[str]]) -> Dict[str, Accession]:
        """ Получаем unused, суммы значений Sc, Precursor Signal и сумму длинн
        последовательностей и подсчитываем количество строк с одинаковым
        Accession для каждого Accession """

        accessions: Dict[str, Accession] = {}
        i = 0
        while i < len(peptideTable[self.columnNames.accession]):
            curAccession = (
                peptideTable[self.columnNames.accession][i].split(sep=';')[0])
            if curAccession not in accessions:
                accessions[curAccession] = Accession(name=curAccession)
            accessions[curAccession].Counts += 1
            accessions[curAccession].Unused = Decimal(
                peptideTable[self.columnNames.unused][i])
            accessions[curAccession].ScSumm += Decimal(
                peptideTable[self.columnNames.sc][i])
            accessions[curAccession].PSignalSumm += (
                Decimal(peptideTable[self.columnNames.precursorSignal][i]) if
                peptideTable[self.columnNames.precursorSignal][i] != ''
                else Decimal(0))
            accessions[curAccession].SeqlenSumm += (
                len(peptideTable[self.columnNames.sequence][i]))
            i += 1
        return accessions

    def CalculateNormParamsForAccessions(
            self,
            accessions: Dict[str, Accession],
            seqences: Dict[str, Sequence]):

        for accession, curAccession in accessions.items():
            curAccession.ScNorm = curAccession.ScSumm / seqences[accession].len
            curAccession.PSignalNorm = (
                curAccession.PSignalSumm / seqences[accession].len)

    def GenerateAccessionsBunchOverAllTables(
            self) -> Dict[str, Dict[str, Accession]]:

        accessions: Dict[str, Dict[str, Accession]] = {}
        for tableName, table in self.accessionsPerTable.items():
            for accessionName, accession in table.items():
                if accessionName not in accessions:
                    accessions[accessionName] = {}
                accessions[accessionName][tableName] = accession
        return accessions

    def RemoveAccessionFromAllTables(self, accession: str) -> None:

        for table in self.accessionsPerTable.values():
            table.pop(accession, None)

    def SetColumnNames(self, columnNames: ColumnNames):

        self.columnNames = columnNames
