from typing import Dict, List
from Classes.Accession import Accession
from Classes.Sequence import Sequence
from Classes.PeptideTables import PeptideTables


class AccessionTables:

    accessionsPerTable: Dict[str, Dict[str, Accession]]
    sortedTableNums: List[str]

    def __init__(self,
                 seqDB: Dict[str, Sequence],
                 peptideTables: PeptideTables) -> None:
        self.proteinReplacements = None
        self.GetAccessionsPerTable(seqDB, peptideTables)

    def GetAccessionsFromTable(
            self,
            peptideTable: Dict[str, List[str]]) -> Dict[str, Accession]:
        """ Получаем unused, суммы значений Sc, Precursor Signal и сумму длинн
        последовательностей и подсчитываем количество строк с одинаковым
        Accession для каждого Accession """

        accessions: Dict[str, Accession] = {}
        i = 0
        while i < len(peptideTable["Accessions"]):
            curAccession = peptideTable["Accessions"][i].split(sep=';')[0]
            if curAccession not in accessions:
                accessions[curAccession] = Accession(name=curAccession)
            accessions[curAccession].Counts += 1
            accessions[curAccession].Unused = float(peptideTable["Unused"][i])
            accessions[curAccession].ScSumm += float(peptideTable["Sc"][i])
            accessions[curAccession].PSignalSumm += float(
                peptideTable["PrecursorSignal"][i])
            accessions[curAccession].SeqlenSumm += (
                len(peptideTable["Sequence"][i]))
            i += 1
        return accessions

    def CalculateNormParamsForAccessions(self: "AccessionTables",
                                         accessions: Dict[str, Accession],
                                         seqences: Dict[str, Sequence]):
        for accession in accessions:
            curAccession = accessions[accession]
            curAccession.ScNorm = curAccession.ScSumm / seqences[accession].len
            curAccession.PSignalNorm = (
                curAccession.PSignalSumm / seqences[accession].len)

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
